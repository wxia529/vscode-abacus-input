/**
 * Completion Provider for ABACUS INPUT files
 * Provides intelligent keyword completion with enum value support
 */

import * as vscode from 'vscode';
import { AbacusDataProvider, CompletionItem } from './dataProvider';

export class AbacusCompletionProvider implements vscode.CompletionItemProvider {
  private dataProvider: AbacusDataProvider;

  constructor(dataProvider: AbacusDataProvider) {
    this.dataProvider = dataProvider;
  }

  public provideCompletionItems(
    document: vscode.TextDocument,
    position: vscode.Position,
    _token: vscode.CancellationToken,
    context: vscode.CompletionContext
  ): vscode.CompletionItem[] | vscode.CompletionList {
    const config = vscode.workspace.getConfiguration('abacusInput');
    if (!config.get('completion.enabled', true)) {
      return [];
    }

    const linePrefix = document.lineAt(position).text.slice(0, position.character);
    
    // Check if we're completing a value (keyword already typed with space and optional partial value)
    // Match patterns like: "relax_new " or "relax_new T" or "  calculation scf"
    const valueCompletionMatch = linePrefix.match(/^\s*(\w+)\s+(\w*)$/);
    const isValueCompletion = !!valueCompletionMatch;
    
    // Get the keyword and partial value if we're in value completion mode
    let currentKeyword: string | undefined;
    let partialValue: string = '';
    if (isValueCompletion) {
      currentKeyword = valueCompletionMatch![1];
      partialValue = valueCompletionMatch![2] || '';
    }
    
    // Only trigger on letter input or explicit trigger (Ctrl+Space)
    if (!context.triggerCharacter && context.triggerKind !== vscode.CompletionTriggerKind.Invoke) {
      // Check if we're at the start of a line or after whitespace (keyword position)
      if (!/^\s*[a-zA-Z_]*$/.test(linePrefix) && !/^\s*$/.test(linePrefix)) {
        return [];
      }
    }

    const completionItems = this.dataProvider.getCompletionItems();
    const items: vscode.CompletionItem[] = [];

    // If in value completion mode, only show value suggestions for the keyword
    if (isValueCompletion && currentKeyword) {
      const keywordItem = completionItems.find(item => item.label === currentKeyword);
      if (keywordItem && keywordItem.enumValues && keywordItem.enumValues.length > 0) {
        // Filter enum values by partial value
        const filteredValues = keywordItem.enumValues.filter(v => 
          v.toLowerCase().startsWith(partialValue.toLowerCase())
        );
        
        // Create value completion items
        for (const value of filteredValues) {
          const valueItem = new vscode.CompletionItem(
            value,
            vscode.CompletionItemKind.Value
          );
          valueItem.detail = keywordItem.detail;
          valueItem.documentation = new vscode.MarkdownString(
            `Value for \`${keywordItem.label}\`\n\n${keywordItem.documentation}`
          );
          // Insert the full value, replacing the partial value
          valueItem.insertText = value;
          valueItem.range = new vscode.Range(
            new vscode.Position(position.line, position.character - partialValue.length),
            position
          );
          valueItem.filterText = value;
          valueItem.sortText = keywordItem.section + '_' + keywordItem.label + '_' + value;
          items.push(valueItem);
        }
        return items;
      }
      // If no enum values, return empty (user needs to type the value)
      return [];
    }

    // Group items by section for better organization
    const sectionGroups = new Map<string, vscode.CompletionItem[]>();

    for (const item of completionItems) {
      const completionItem = this.createCompletionItem(item);
      
      if (!sectionGroups.has(item.section)) {
        sectionGroups.set(item.section, []);
      }
      sectionGroups.get(item.section)!.push(completionItem);
    }

    // Sort sections and add to result
    const sortedSections = Array.from(sectionGroups.keys()).sort();
    for (const section of sortedSections) {
      const groupItems = sectionGroups.get(section)!;
      // Sort items within section by label
      groupItems.sort((a, b) => String(a.label).localeCompare(String(b.label)));
      items.push(...groupItems);
    }

    return items;
  }

  private createCompletionItem(item: CompletionItem): vscode.CompletionItem {
    const completionItem = new vscode.CompletionItem(
      item.label,
      this.getCompletionItemKind(item.kind)
    );

    completionItem.detail = item.detail;
    completionItem.documentation = new vscode.MarkdownString(item.documentation);
    completionItem.filterText = item.label;
    completionItem.sortText = item.section + '_' + item.label;

    // Insert keyword with trailing space, then user types or selects value
    completionItem.insertText = new vscode.SnippetString(`${item.label} `);

    return completionItem;
  }

  private getCompletionItemKind(kind: string): vscode.CompletionItemKind {
    switch (kind.toLowerCase()) {
      case 'keyword':
        return vscode.CompletionItemKind.Keyword;
      case 'function':
        return vscode.CompletionItemKind.Function;
      case 'variable':
        return vscode.CompletionItemKind.Variable;
      case 'constant':
        return vscode.CompletionItemKind.Constant;
      default:
        return vscode.CompletionItemKind.Text;
    }
  }

  public resolveCompletionItem?(
    item: vscode.CompletionItem,
    _token: vscode.CancellationToken
  ): vscode.CompletionItem {
    // Additional resolution if needed
    return item;
  }
}
