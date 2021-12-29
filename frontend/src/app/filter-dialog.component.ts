import { Component, OnInit } from "@angular/core";
import { ApiService } from "./api.service";
import {
  AbstractControl,
  FormControl,
  FormGroup,
  ValidationErrors,
  Validators,
} from "@angular/forms";
import { map, Observable } from "rxjs";
import { MatDialogRef } from "@angular/material/dialog";

@Component({
  selector: "app-filter-dialog",
  templateUrl: "filter-dialog.component.html",
  styleUrls: ["filter-dialog.component.scss"],
  // changeDetection: ChangeDetectionStrategy.OnPush,
})
export class FilterDialogComponent implements OnInit {
  public filterForm: FormGroup = new FormGroup({
    smarts: new FormControl("", {
      validators: [Validators.required],
      asyncValidators: [(control) => this.validateSubstructure(control)],
      updateOn: "blur",
    }),
  });

  constructor(
    private apiService: ApiService,
    private dialogRef: MatDialogRef<FilterDialogComponent>
  ) {}

  ngOnInit(): void {}

  validateSubstructure(
    control: AbstractControl
  ): Observable<ValidationErrors | null> {
    console.log(this, control);

    return this.apiService.isSubstructureValid(control.value).pipe(
      map((isValid) => {
        return isValid ? null : { invalid: true };
      })
    );
  }

  getSubstructureErrorMessage() {
    if (this.filterForm.get("smarts")?.hasError("required"))
      return "enter a value";
    return this.filterForm.get("smarts")?.hasError("invalid")
      ? "invalid SMARTS"
      : "";
  }

  onSubmit() {
    if (this.filterForm.invalid) return;

    console.log(this.filterForm.value);
    this.dialogRef.close(this.filterForm.value);
  }
}
