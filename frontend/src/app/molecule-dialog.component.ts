import {
  ChangeDetectionStrategy,
  ChangeDetectorRef,
  Component,
  Inject,
  OnDestroy,
  OnInit,
} from "@angular/core";
import { ApiService } from "./api.service";
import { MAT_DIALOG_DATA } from "@angular/material/dialog";
import { BASE_API_URL } from "./app.module";
import {
  GETMoleculeResponse,
  GETMoleculeResponseBase,
  MoleculeDescriptors,
} from "./api.interface";
import { MatTooltip } from "@angular/material/tooltip";
import { Subscription } from "rxjs";

@Component({
  selector: "app-filter-dialog",
  templateUrl: "molecule-dialog.component.html",
  styleUrls: ["molecule-dialog.component.scss"],
  changeDetection: ChangeDetectionStrategy.OnPush,
})
export class MoleculeDialogComponent implements OnInit, OnDestroy {
  public readonly BASE_API_URL: string = BASE_API_URL;
  public descriptors?: MoleculeDescriptors = undefined;
  private moleculeSubscription$: Subscription;

  constructor(
    private apiService: ApiService,
    private changeRef: ChangeDetectorRef,
    @Inject(MAT_DIALOG_DATA) public data: { molecule: GETMoleculeResponseBase }
  ) {
    console.log(this.data.molecule.self);

    this.moleculeSubscription$ = this.apiService
      .getEndpoint<GETMoleculeResponse>(this.data.molecule.self)
      .subscribe((response) => {
        this.descriptors = response.descriptors;
        this.changeRef.detectChanges();
      });
  }

  ngOnInit(): void {}

  ngOnDestroy() {
    this.moleculeSubscription$.unsubscribe();
  }

  onCopied(tooltip: MatTooltip, success: boolean) {
    tooltip.disabled = false;
    tooltip.message = success ? "Copied!" : "Error :(";
    tooltip.show();
    setTimeout(() => {
      tooltip.hide();
      tooltip.disabled = true;
    }, 500);
  }
}
