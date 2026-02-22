/**
 * Hover Provider for ABACUS INPUT files
 * Displays documentation when hovering over keywords
 */
import * as vscode from 'vscode';
import { AbacusDataProvider } from './dataProvider';
export declare class AbacusHoverProvider implements vscode.HoverProvider {
    private dataProvider;
    constructor(dataProvider: AbacusDataProvider);
    provideHover(document: vscode.TextDocument, position: vscode.Position, _token: vscode.CancellationToken): vscode.Hover | undefined;
    private getLineDiagnostics;
    private buildMarkdownContent;
    private createKeywordHover;
    private createConstraintHover;
    private createValueHover;
    private formatDiagnostics;
    private createDiagnosticsHover;
}
//# sourceMappingURL=hoverProvider.d.ts.map