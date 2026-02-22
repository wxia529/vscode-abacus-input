/**
 * Diagnostic Provider for ABACUS INPUT files
 * Validates parameters against constraints and dependency rules
 */
import * as vscode from 'vscode';
import { AbacusDataProvider } from './dataProvider';
export declare class AbacusDiagnosticProvider {
    private dataProvider;
    private diagnosticCollection;
    constructor(dataProvider: AbacusDataProvider);
    dispose(): void;
    diagnose(document: vscode.TextDocument, token: vscode.CancellationToken): Promise<void>;
    private parseDocument;
    private parseLine;
    private validateLine;
    private validateType;
    private validateEnum;
    private validateRange;
    private validateDependencies;
    private validateConditionalRules;
    private evaluateCondition;
    private valuesMatch;
    private getDiagnosticSeverity;
    private createDiagnostic;
}
//# sourceMappingURL=diagnosticProvider.d.ts.map