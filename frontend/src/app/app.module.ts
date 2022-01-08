import { NgModule } from "@angular/core";
import { BrowserModule } from "@angular/platform-browser";
import { BrowserAnimationsModule } from "@angular/platform-browser/animations";
import { HttpClientModule } from "@angular/common/http";
import { FormsModule, ReactiveFormsModule } from "@angular/forms";

import { MaterialModule } from "./material.module";

import { AppRoutingModule } from "./app-routing.module";
import { AppComponent } from "./app.component";
import { PaginatorComponent } from "./paginator.component";
import { FilterDialogComponent } from "./filter-dialog.component";
import { MoleculeDialogComponent } from "./molecule-dialog.component";
import { ClipboardModule } from "@angular/cdk/clipboard";

declare var baseAPIUrl: string;
export const BASE_API_URL = baseAPIUrl;

@NgModule({
  declarations: [
    AppComponent,
    PaginatorComponent,
    FilterDialogComponent,
    MoleculeDialogComponent,
  ],
  imports: [
    MaterialModule,

    BrowserModule,
    BrowserAnimationsModule,

    HttpClientModule,

    ClipboardModule,

    FormsModule,
    ReactiveFormsModule,

    AppRoutingModule,
  ],
  providers: [],
  bootstrap: [AppComponent],
})
export class AppModule {}
