"use strict";
/**
 * Completion Provider for ABACUS INPUT files
 * Provides intelligent keyword completion with enum value support
 */
var __createBinding = (this && this.__createBinding) || (Object.create ? (function(o, m, k, k2) {
    if (k2 === undefined) k2 = k;
    var desc = Object.getOwnPropertyDescriptor(m, k);
    if (!desc || ("get" in desc ? !m.__esModule : desc.writable || desc.configurable)) {
      desc = { enumerable: true, get: function() { return m[k]; } };
    }
    Object.defineProperty(o, k2, desc);
}) : (function(o, m, k, k2) {
    if (k2 === undefined) k2 = k;
    o[k2] = m[k];
}));
var __setModuleDefault = (this && this.__setModuleDefault) || (Object.create ? (function(o, v) {
    Object.defineProperty(o, "default", { enumerable: true, value: v });
}) : function(o, v) {
    o["default"] = v;
});
var __importStar = (this && this.__importStar) || (function () {
    var ownKeys = function(o) {
        ownKeys = Object.getOwnPropertyNames || function (o) {
            var ar = [];
            for (var k in o) if (Object.prototype.hasOwnProperty.call(o, k)) ar[ar.length] = k;
            return ar;
        };
        return ownKeys(o);
    };
    return function (mod) {
        if (mod && mod.__esModule) return mod;
        var result = {};
        if (mod != null) for (var k = ownKeys(mod), i = 0; i < k.length; i++) if (k[i] !== "default") __createBinding(result, mod, k[i]);
        __setModuleDefault(result, mod);
        return result;
    };
})();
Object.defineProperty(exports, "__esModule", { value: true });
exports.AbacusCompletionProvider = void 0;
const vscode = __importStar(require("vscode"));
class AbacusCompletionProvider {
    constructor(dataProvider) {
        this.dataProvider = dataProvider;
    }
    provideCompletionItems(document, position, _token, context) {
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
        let currentKeyword;
        let partialValue = '';
        if (isValueCompletion) {
            currentKeyword = valueCompletionMatch[1];
            partialValue = valueCompletionMatch[2] || '';
        }
        // Only trigger on letter input or explicit trigger (Ctrl+Space)
        if (!context.triggerCharacter && context.triggerKind !== vscode.CompletionTriggerKind.Invoke) {
            // Check if we're at the start of a line or after whitespace (keyword position)
            if (!/^\s*[a-zA-Z_]*$/.test(linePrefix) && !/^\s*$/.test(linePrefix)) {
                return [];
            }
        }
        const completionItems = this.dataProvider.getCompletionItems();
        const items = [];
        // If in value completion mode, only show value suggestions for the keyword
        if (isValueCompletion && currentKeyword) {
            const keywordItem = completionItems.find(item => item.label === currentKeyword);
            if (keywordItem && keywordItem.enumValues && keywordItem.enumValues.length > 0) {
                // Filter enum values by partial value
                const filteredValues = keywordItem.enumValues.filter(v => v.toLowerCase().startsWith(partialValue.toLowerCase()));
                // Create value completion items
                for (const value of filteredValues) {
                    const valueItem = new vscode.CompletionItem(value, vscode.CompletionItemKind.Value);
                    valueItem.detail = keywordItem.detail;
                    valueItem.documentation = new vscode.MarkdownString(`Value for \`${keywordItem.label}\`\n\n${keywordItem.documentation}`);
                    // Insert the full value, replacing the partial value
                    valueItem.insertText = value;
                    valueItem.range = new vscode.Range(new vscode.Position(position.line, position.character - partialValue.length), position);
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
        const sectionGroups = new Map();
        for (const item of completionItems) {
            const completionItem = this.createCompletionItem(item);
            if (!sectionGroups.has(item.section)) {
                sectionGroups.set(item.section, []);
            }
            sectionGroups.get(item.section).push(completionItem);
        }
        // Sort sections and add to result
        const sortedSections = Array.from(sectionGroups.keys()).sort();
        for (const section of sortedSections) {
            const groupItems = sectionGroups.get(section);
            // Sort items within section by label
            groupItems.sort((a, b) => String(a.label).localeCompare(String(b.label)));
            items.push(...groupItems);
        }
        return items;
    }
    createCompletionItem(item) {
        const completionItem = new vscode.CompletionItem(item.label, this.getCompletionItemKind(item.kind));
        completionItem.detail = item.detail;
        completionItem.documentation = new vscode.MarkdownString(item.documentation);
        completionItem.filterText = item.label;
        completionItem.sortText = item.section + '_' + item.label;
        // Insert keyword with trailing space, then user types or selects value
        completionItem.insertText = new vscode.SnippetString(`${item.label} `);
        return completionItem;
    }
    getCompletionItemKind(kind) {
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
    resolveCompletionItem(item, _token) {
        // Additional resolution if needed
        return item;
    }
}
exports.AbacusCompletionProvider = AbacusCompletionProvider;
//# sourceMappingURL=completionProvider.js.map