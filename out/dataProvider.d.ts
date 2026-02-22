/**
 * Data loader for ABACUS INPUT JSON files
 */
import * as vscode from 'vscode';
export interface CompletionItem {
    label: string;
    detail: string;
    documentation: string;
    section: string;
    insertText: string;
    insertTextFormat: number;
    kind: string;
    enumValues: string[];
}
export interface Constraint {
    type: string;
    enum?: string[];
    default?: string;
    unit: string;
    description: string;
    exclusiveMinimum?: number;
    exclusiveMaximum?: number;
    minimum?: number;
    maximum?: number;
    availability?: string;
    mdType?: string;
    section?: string;
}
export interface DiagnosticRule {
    keyword: string;
    severity: string;
    message: string;
    condition?: {
        keyword: string;
        operator: string;
        value: string | number | boolean | string[];
    };
}
export interface Snippet {
    prefix: string;
    body: string;
    description: string;
}
export interface Dependency {
    requires: Record<string, string | number | boolean>;
}
export declare class AbacusDataProvider {
    private completionData;
    private constraintsData;
    private diagnosticsData;
    private snippetsData;
    private dependenciesData;
    private dataPath;
    constructor(context: vscode.ExtensionContext);
    loadData(): Promise<void>;
    private loadJsonFile;
    getCompletionItems(): CompletionItem[];
    getConstraint(keyword: string): Constraint | undefined;
    getDiagnosticRules(): DiagnosticRule[];
    getSnippets(): Record<string, Snippet>;
    getDependency(keyword: string): Dependency | undefined;
    getAllKeywords(): string[];
    getCompletionByLabel(label: string): CompletionItem | undefined;
}
//# sourceMappingURL=dataProvider.d.ts.map