import functools
import os
import sqlite3
from typing import Any, Iterable, List, NamedTuple, Optional, Tuple, Union

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from splore.models import Cursor, Page, RangeFilter, SMARTSFilter, SortBy

_REVERSE_SORT_BY = {None: None, "asc": "desc", "desc": "asc"}

PER_PAGE_DEFAULT = 200


def _sort_by_to_sql(sort_by: List[SortBy], reverse: bool) -> str:

    statements = [
        " ".join((column, direction if not reverse else _REVERSE_SORT_BY[direction]))
        for column, direction in sort_by
    ]

    return f"order by {','.join(statements)}" if len(statements) > 0 else ""


def _filters_to_sql(filters: List[Union[RangeFilter, SMARTSFilter]]) -> str:

    statements = []

    for column_filter in filters:

        if isinstance(column_filter, RangeFilter):

            if column_filter.lt is not None:
                statements.append(f"{column_filter.column} < {column_filter.lt}")
            if column_filter.le is not None:
                statements.append(f"{column_filter.column} <= {column_filter.le}")
            if column_filter.gt is not None:
                statements.append(f"{column_filter.column} > {column_filter.gt}")
            if column_filter.ge is not None:
                statements.append(f"{column_filter.column} >= {column_filter.ge}")

        elif isinstance(column_filter, SMARTSFilter):
            statements.append(f"smarts_match(smiles,'{column_filter.smarts}')")

        else:
            raise NotImplementedError()

    return " and ".join(statements)


@functools.lru_cache(16384)
def _filter_by_pattern(smiles: str, pattern: str) -> bool:

    molecule: Chem.Mol = Chem.MolFromSmiles(smiles)
    q_mol: Chem.Mol = Chem.MolFromSmarts(pattern)

    return molecule.HasSubstructMatch(q_mol)


class SploreDBRow(NamedTuple):

    smiles: str

    weight: float
    n_heavy_atoms: int

    n_aliphatic_carbocycles: int
    n_aliphatic_heterocycles: int

    n_aromatic_carbocycles: int
    n_aromatic_heterocycles: int

    n_rotatable_bonds: int

    n_h_bond_acceptors: int
    n_h_bond_donors: int

    topological_polar_surface_area: float


class SploreDBPage:
    @property
    def has_next(self) -> bool:
        return len(self.rows) > 0 and self._cursor_next is not None

    @property
    def has_prev(self) -> bool:
        return len(self.rows) > 0 and self._cursor_prev is not None

    @property
    def next(self) -> Page:
        return (self._cursor_end or self._cursor_prev), "next"

    @property
    def prev(self) -> Page:
        return (self._cursor_start or self._cursor_next), "prev"

    @property
    def current(self) -> Page:
        if self.backwards:
            return self._cursor_next, "prev"
        else:
            return self._cursor_prev, "next"

    def __init__(
        self,
        cursor: Cursor,
        backwards: bool,
        per_page: int,
        keys: List[Cursor],
        rows: List[Tuple[Any]],
    ):

        self.per_page = per_page
        self.backwards = backwards

        # try and retrieve an extra row to see if we can go further forward / back
        self.rows = rows[:per_page]

        extra_keys = keys[per_page:]
        keys = keys[:per_page]

        cursors = (
            # current =
            cursor,
            # start and end
            None if len(keys) == 0 else keys[0],
            None if len(keys) == 0 else keys[-1],
            # next
            None if len(extra_keys) == 0 else extra_keys[0],
        )

        if backwards:  # make sure to flip the results if getting previous
            self.rows = self.rows[::-1]
            cursors = cursors[::-1]

        (
            self._cursor_prev,
            self._cursor_start,
            self._cursor_end,
            self._cursor_next,
        ) = cursors

    def __len__(self):
        return len(self.rows)


class SploreDB:
    """A wrapper around a SQLite database that can store and be queried for RDKit
    molecules
    """

    @property
    def n_molecules(self) -> int:

        if self._n_molecules is None:

            self._n_molecules = self._connection.execute(
                "select count(*) from molecules"
            ).fetchone()[0]

        return self._n_molecules

    def __init__(self, file_path: Union[os.PathLike, str], clear_existing: bool = True):

        self._file_path = file_path

        self._connection = sqlite3.connect(file_path)
        self._connection.create_function(
            "smarts_match", 2, _filter_by_pattern, deterministic=True
        )

        self._create_schema()

        if clear_existing:
            self.clear()

        self._n_molecules = None

    def __del__(self):

        if self._connection is not None:
            self._connection.close()

    def _create_schema(self):

        with self._connection:

            self._connection.execute(
                "create table if not exists info (version integer)"
            )
            db_info = self._connection.execute("select * from info").fetchall()

            if len(db_info) == 0:
                self._connection.execute("insert into info values(1)")
            else:
                assert len(db_info) == 1 and db_info[0] == (1,)

            python_to_sql_type = {str: "text", int: "integer", float: "real"}

            columns = ", ".join(
                f"{field} {python_to_sql_type[field_type]}"
                for field, field_type in SploreDBRow.__annotations__.items()
            )

            self._connection.execute(
                f"create table if not exists molecules ( {columns} )"
            )

            for column in SploreDBRow.__annotations__:

                if column in {"smiles"}:
                    continue

                self._connection.execute(
                    f"create index if not exists ix_{column} on molecules({column})"
                )

            self._connection.execute("pragma optimize")

    def create(self, molecules: Iterable[Chem.Mol]):

        with self._connection:

            attributes = ", ".join(attr for attr in SploreDBRow.__annotations__)
            value_stub = ", ".join("?" * len(SploreDBRow.__annotations__))

            self._connection.executemany(
                f"insert into molecules ({attributes}) values ({value_stub})",
                (
                    SploreDBRow(
                        Chem.MolToSmiles(molecule),
                        rdMolDescriptors.CalcExactMolWt(molecule),
                        rdMolDescriptors.CalcNumHeavyAtoms(molecule),
                        rdMolDescriptors.CalcNumAliphaticCarbocycles(molecule),
                        rdMolDescriptors.CalcNumAliphaticHeterocycles(molecule),
                        rdMolDescriptors.CalcNumAromaticCarbocycles(molecule),
                        rdMolDescriptors.CalcNumAromaticHeterocycles(molecule),
                        rdMolDescriptors.CalcNumRotatableBonds(molecule),
                        rdMolDescriptors.CalcNumHBA(molecule),
                        rdMolDescriptors.CalcNumHBD(molecule),
                        rdMolDescriptors.CalcTPSA(molecule, includeSandP=True),
                    )
                    for molecule in molecules
                ),
            )
            self._connection.execute("pragma optimize")

        self._n_molecules = None

    def read_all(
        self,
        page: Page = (None, "next"),
        per_page: int = PER_PAGE_DEFAULT,
        sort_by: Optional[SortBy] = None,
        filters: Optional[List[Union[RangeFilter, SMARTSFilter]]] = None,
    ) -> SploreDBPage:

        cursor, move_to = page
        backwards = {"prev": True, "next": False}[move_to]

        order_by_terms = ([] if not sort_by else [sort_by]) + [("ROWID", "asc")]
        order_by_sql = _sort_by_to_sql(order_by_terms, reverse=backwards)

        filters = filters if filters else []
        filters_sql = _filters_to_sql(filters)

        if cursor is not None:

            iterable_cursor = (cursor,) if not isinstance(cursor, tuple) else cursor

            lhs = tuple(
                column
                if (
                    (not backwards and order_by == "asc")
                    or (backwards and order_by == "desc")
                )
                else value
                for value, (column, order_by) in zip(iterable_cursor, order_by_terms)
            )
            rhs = tuple(
                column
                if not (
                    (not backwards and order_by == "asc")
                    or (backwards and order_by == "desc")
                )
                else value
                for value, (column, order_by) in zip(iterable_cursor, order_by_terms)
            )
            assert len(lhs) == len(rhs) and len(lhs) == len(iterable_cursor)

            if len(lhs) == 1:
                page_sql = f"{lhs[0]} > {rhs[0]}"
            else:
                page_sql = f"({','.join(map(str, lhs))}) > ({','.join(map(str, rhs))})"

            filters_sql = (
                f"{page_sql} and {filters_sql}" if filters_sql != "" else page_sql
            )

        where_sql = f"where {filters_sql}" if len(filters_sql) > 0 else ""

        required_fields = [column for column, _ in order_by_terms] + ["smiles"]

        for range_filter in filters:
            if (
                not isinstance(filter, RangeFilter)
                or range_filter.column in required_fields
            ):
                continue
            required_fields.append(range_filter.column)

        select_sql = ",".join(required_fields)

        limit = per_page + 1

        statement = (
            f"select {select_sql} from molecules "
            f"{where_sql} "
            f"{order_by_sql} "
            f"limit {limit} "
        )

        rows = self._connection.execute(statement).fetchall()
        keys = [row[: len(order_by_terms)] for row in rows]

        return SploreDBPage(cursor, backwards, per_page, keys, rows)

    def read(self, molecule_id: int) -> SploreDBRow:

        return SploreDBRow(
            *self._connection.execute(
                f"select * from molecules where ROWID = {molecule_id}"
            ).fetchone()
        )

    def clear(self):

        with self._connection:
            self._connection.execute("delete from molecules")

        self._n_molecules = None
