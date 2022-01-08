import {
  ChangeDetectionStrategy,
  Component,
  ElementRef,
  Inject,
  OnInit,
} from "@angular/core";
import { ApiService } from "./api.service";
import { MAT_DIALOG_DATA } from "@angular/material/dialog";
import { BASE_API_URL } from "./app.module";
import { GETMoleculeResponse } from "./api.interface";
import { MatTooltip } from "@angular/material/tooltip";

@Component({
  selector: "app-filter-dialog",
  templateUrl: "molecule-dialog.component.html",
  styleUrls: ["molecule-dialog.component.scss"],
  changeDetection: ChangeDetectionStrategy.OnPush,
})
export class MoleculeDialogComponent implements OnInit {
  public readonly BASE_API_URL: string = BASE_API_URL;

  constructor(
    private apiService: ApiService,
    private element: ElementRef,
    @Inject(MAT_DIALOG_DATA) public data: { molecule: GETMoleculeResponse }
  ) {}

  ngOnInit(): void {}

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
