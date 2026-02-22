/**
 * Completion Provider for ABACUS INPUT files
 * Provides intelligent keyword completion with enum value support
 */
import * as vscode from 'vscode';
import { AbacusDataProvider } from './dataProvider';
export declare class AbacusCompletionProvider implements vscode.CompletionItemProvider {
    private dataProvider;
    constructor(dataProvider: AbacusDataProvider);
    provideCompletionItems(document: vscode.TextDocument, position: vscode.Position, _token: vscode.CancellationToken, context: vscode.CompletionContext): vscode.CompletionItem[] | vscode.CompletionList;
    private createCompletionItem;
    private getCompletionItemKind;
    resolveCompletionItem?(item: vscode.CompletionItem, _token: vscode.CancellationToken): vscode.CompletionItem;
}
//# sourceMappingURL=completionProvider.d.ts.map