"use strict";
/**
 * Hover Provider for ABACUS INPUT files
 * Displays documentation when hovering over keywords
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
exports.AbacusHoverProvider = void 0;
const vscode = __importStar(require("vscode"));
class AbacusHoverProvider {
    constructor(dataProvider) {
        this.dataProvider = dataProvider;
    }
    provideHover(document, position, _token) {
        const config = vscode.workspace.getConfiguration('abacusInput');
        if (!config.get('hover.enabled', true)) {
            return undefined;
        }
        const wordRange = document.getWordRangeAtPosition(position);
        if (!wordRange) {
            return undefined;
        }
        const word = document.getText(wordRange);
        const lineText = document.lineAt(position).text;
        const lineDiagnostics = this.getLineDiagnostics(document, position);
        // Check if this is a value hover (previous word is a keyword)
        const linePrefix = lineText.slice(0, wordRange.start.character);
        const keywordMatch = linePrefix.match(/\s*(\w+)\s*$/);
        if (keywordMatch && !/^\s*$/.test(linePrefix)) {
            const prevKeyword = keywordMatch[1];
            const constraint = this.dataProvider.getConstraint(prevKeyword);
            if (constraint) {
                return this.createValueHover(constraint, word, lineDiagnostics);
            }
            return undefined;
        }
        // Check for keyword completion item
        const completionItem = this.dataProvider.getCompletionByLabel(word);
        if (completionItem) {
            return this.createKeywordHover(completionItem, lineDiagnostics);
        }
        // Check for constraint
        const constraint = this.dataProvider.getConstraint(word);
        if (constraint) {
            return this.createConstraintHover(word, constraint, lineDiagnostics);
        }
        // Show diagnostics even if no keyword documentation found
        if (lineDiagnostics.length > 0) {
            return this.createDiagnosticsHover(lineDiagnostics);
        }
        return undefined;
    }
    getLineDiagnostics(document, position) {
        const line = document.lineAt(position.line);
        const allDiagnostics = vscode.languages.getDiagnostics(document.uri);
        return allDiagnostics.filter(d => d.range.start.line === line.lineNumber);
    }
    buildMarkdownContent(title, type, description, enumValues, defaultValue, unit, availability, section) {
        const markdown = new vscode.MarkdownString();
        markdown.isTrusted = true;
        markdown.supportHtml = true;
        markdown.appendMarkdown(`#### **${title}**\n\n`);
        markdown.appendMarkdown(`*Type:* ${type}\n\n`);
        markdown.appendMarkdown(`${description}\n\n`);
        if (enumValues && enumValues.length > 0) {
            markdown.appendMarkdown('---\n\n');
            markdown.appendMarkdown('**Allowed values:**\n\n');
            markdown.appendMarkdown(enumValues.map(v => `- \`${v}\``).join('  \n'));
            markdown.appendMarkdown('\n\n');
        }
        if (defaultValue !== undefined) {
            markdown.appendMarkdown('---\n\n');
            markdown.appendMarkdown(`**Default:** \`${defaultValue}\`\n\n`);
        }
        if (unit) {
            markdown.appendMarkdown(`**Unit:** ${unit}\n\n`);
        }
        if (availability) {
            markdown.appendMarkdown('---\n\n');
            markdown.appendMarkdown(`**Availability:** ${availability}\n\n`);
        }
        if (section) {
            markdown.appendMarkdown(`---\n\n*Section:* ${section}`);
        }
        return markdown;
    }
    createKeywordHover(item, lineDiagnostics) {
        const markdown = this.buildMarkdownContent(item.label, item.detail, item.documentation, item.enumValues, undefined, undefined, undefined, item.section);
        if (lineDiagnostics.length > 0) {
            const diagnosticContent = this.formatDiagnostics(lineDiagnostics);
            const fullContent = new vscode.MarkdownString();
            fullContent.isTrusted = true;
            fullContent.appendMarkdown(diagnosticContent);
            fullContent.appendMarkdown('\n\n---\n\n');
            fullContent.appendMarkdown(markdown.value);
            return new vscode.Hover(fullContent);
        }
        return new vscode.Hover(markdown);
    }
    createConstraintHover(keyword, constraint, lineDiagnostics) {
        const markdown = this.buildMarkdownContent(keyword, constraint.type, constraint.description, constraint.enum, constraint.default, constraint.unit, constraint.availability);
        if (lineDiagnostics.length > 0) {
            const diagnosticContent = this.formatDiagnostics(lineDiagnostics);
            const fullContent = new vscode.MarkdownString();
            fullContent.isTrusted = true;
            fullContent.appendMarkdown(diagnosticContent);
            fullContent.appendMarkdown('\n\n---\n\n');
            fullContent.appendMarkdown(markdown.value);
            return new vscode.Hover(fullContent);
        }
        return new vscode.Hover(markdown);
    }
    createValueHover(constraint, value, lineDiagnostics) {
        if (!constraint.enum) {
            return undefined;
        }
        const markdown = new vscode.MarkdownString();
        const isValid = constraint.enum.includes(value);
        if (lineDiagnostics.length > 0) {
            markdown.appendMarkdown(this.formatDiagnostics(lineDiagnostics));
            markdown.appendMarkdown('\n\n---\n\n');
        }
        if (isValid) {
            markdown.appendMarkdown(`\`${value}\` is a valid value for **${constraint.type}**\n\n`);
            markdown.appendMarkdown('**All allowed values:**\n\n');
            markdown.appendMarkdown(constraint.enum.map(v => `- \`${v}\``).join('  \n'));
        }
        else {
            markdown.appendMarkdown(`⚠️ \`${value}\` is not a valid value\n\n`);
            markdown.appendMarkdown('**Allowed values:**\n\n');
            markdown.appendMarkdown(constraint.enum.map(v => `- \`${v}\``).join('  \n'));
        }
        return new vscode.Hover(markdown);
    }
    formatDiagnostics(diagnostics) {
        return diagnostics.map(d => {
            const icon = d.severity === vscode.DiagnosticSeverity.Error ? '❌' :
                d.severity === vscode.DiagnosticSeverity.Warning ? '⚠️' : 'ℹ️';
            return `${icon} ${d.message}`;
        }).join('\n\n');
    }
    createDiagnosticsHover(diagnostics) {
        return new vscode.Hover(this.formatDiagnostics(diagnostics));
    }
}
exports.AbacusHoverProvider = AbacusHoverProvider;
//# sourceMappingURL=hoverProvider.js.map