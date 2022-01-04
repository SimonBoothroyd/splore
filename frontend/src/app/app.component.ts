import {
  ChangeDetectionStrategy,
  ChangeDetectorRef,
  Component,
  OnDestroy,
  OnInit,
} from "@angular/core";
import { GETMoleculesResponse } from "./api.interface";
import { BehaviorSubject, mergeMap, Subject, takeUntil } from "rxjs";
import { ApiService } from "./api.service";
import { ActivatedRoute, Router } from "@angular/router";
import { MatDialog } from "@angular/material/dialog";
import { FilterDialogComponent } from "./filter-dialog.component";
import { BASE_API_URL } from "./app.module";

@Component({
  selector: "app-root",
  templateUrl: "app.component.html",
  styleUrls: ["app.component.scss"],
  changeDetection: ChangeDetectionStrategy.OnPush,
})
export class AppComponent implements OnInit, OnDestroy {
  public readonly BASE_API_URL: string = BASE_API_URL;
  response: GETMoleculesResponse | undefined;

  destroy$: Subject<boolean> = new Subject<boolean>();
  isLoading$: BehaviorSubject<boolean> = new BehaviorSubject<boolean>(true);

  public constructor(
    private apiService: ApiService,
    private changeRef: ChangeDetectorRef,
    private route: ActivatedRoute,
    private router: Router,
    private dialog: MatDialog
  ) {}

  public ngOnInit() {
    this.route.queryParamMap
      .pipe(
        mergeMap((params) => {
          this.isLoading$.next(true);

          return this.apiService.getMolecules(
            params.get("page") || undefined,
            params.get("per_page") || undefined,
            params.get("sort_by") || undefined,
            params.has("substr")
              ? [
                  {
                    type: "smarts",
                    smarts: atob(params.get("substr") as string),
                  },
                ]
              : []
          );
        })
      )
      .pipe(takeUntil(this.destroy$))
      .subscribe((response) => {
        this.isLoading$.next(false);

        this.response = response;
        this.changeRef.detectChanges();
      });
  }

  public ngOnDestroy() {
    this.destroy$.next(true);
    this.destroy$.unsubscribe();
  }

  public onPageChanged(endpoint: string | undefined) {
    if (!endpoint) return;

    this.router
      .navigateByUrl(endpoint.replace("/molecules", "/"))
      .catch(console.error);
  }

  public onSort(column: string | undefined, order_by: string | undefined) {
    const httpQueryParams = this.apiService.getMoleculesParams(
      undefined,
      this.response?._metadata.per_page,
      column && order_by ? [column, order_by] : undefined,
      this.response?._metadata.filters
    );

    const queryParams: { [key: string]: string } = {};
    httpQueryParams.keys().forEach((key) => {
      queryParams[key] = httpQueryParams.get(key) as string;
    });

    this.router.navigate(["/"], { queryParams }).catch(console.error);
  }

  public onFilterClicked() {
    const dialogRef = this.dialog.open(FilterDialogComponent);
    dialogRef.afterClosed().subscribe((result) => this.onFilter(result));
  }

  public onFilter(result?: any) {
    if (!result) return;

    console.log(this.response);

    const httpQueryParams = this.apiService.getMoleculesParams(
      undefined,
      this.response?._metadata.per_page,
      this.response?._metadata.sort_by,
      result.smarts ? [{ type: "smarts", smarts: result.smarts }] : []
    );

    const queryParams: { [key: string]: string } = {};
    httpQueryParams.keys().forEach((key) => {
      queryParams[key] = httpQueryParams.get(key) as string;
    });

    this.router.navigate(["/"], { queryParams }).catch(console.error);
  }
}
