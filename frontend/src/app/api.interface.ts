export interface RangeFilter {
  type: string;

  column: string;

  lt?: number;
  le?: number;

  gt?: number;
  ge?: number;
}

export interface SMARTSFilter {
  type: string;

  smarts: string;
}

export interface MoleculeDescriptors {
  weight: number;
  n_heavy_atoms: number;

  n_aliphatic_carbocycles: number;
  n_aliphatic_heterocycles: number;

  n_aromatic_carbocycles: number;
  n_aromatic_heterocycles: number;

  n_rotatable_bonds: number;

  n_h_bond_acceptors: number;
  n_h_bond_donors: number;

  topological_polar_surface_area: number;
}

export interface GETMoleculeResponseLinks {
  img: string;
}

export interface GETMoleculeResponseBase {
  self: string;
  id: number;

  smiles: string;

  _links: GETMoleculeResponseLinks;
}

export interface GETMoleculeResponse extends GETMoleculeResponseBase {
  self: string;
  id: number;

  smiles: string;
  descriptors: MoleculeDescriptors;

  _links: GETMoleculeResponseLinks;
}

export interface GETMoleculesResponseMetadata {
  cursor?: string;
  move_to: string;

  per_page: number;

  sort_by?: [string, string];
  filters: (SMARTSFilter | RangeFilter)[];
}

export interface GETMoleculesResponseLinks {
  self: string;

  first: string;
  prev: string;
  next: string;
  last: string;
}

export interface GETMoleculesResponse {
  _metadata: GETMoleculesResponseMetadata;
  _links: GETMoleculesResponseLinks;

  contents: GETMoleculeResponseBase[];
}
