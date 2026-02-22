"use strict";
/**
 * ABACUS INPUT Helper Extension for VSCode
 * Provides intelligent completion, validation, and diagnostics for ABACUS INPUT files
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
exports.activate = activate;
exports.deactivate = deactivate;
const vscode = __importStar(require("vscode"));
const dataProvider_1 = require("./dataProvider");
const completionProvider_1 = require("./completionProvider");
const hoverProvider_1 = require("./hoverProvider");
const diagnosticProvider_1 = require("./diagnosticProvider");
let diagnosticProvider;
let dataLoaded = false;
/**
 * Run diagnostics on all open abacus-input documents
 */
function diagnoseAllOpenDocuments() {
    vscode.workspace.textDocuments.forEach((document) => {
        if (document.languageId === 'abacus-input') {
            const tokenSource = new vscode.CancellationTokenSource();
            diagnosticProvider?.diagnose(document, tokenSource.token);
        }
    });
}
function activate(context) {
    console.log('ABACUS INPUT Helper extension is now active');
    const dataProvider = new dataProvider_1.AbacusDataProvider(context);
    dataProvider.loadData().then(() => {
        dataLoaded = true;
        console.log('ABACUS data loaded successfully');
        diagnoseAllOpenDocuments();
    }).catch((error) => {
        console.error('Failed to load ABACUS data:', error);
        vscode.window.showErrorMessage('Failed to load ABACUS INPUT data. Please check the data path configuration.');
    });
    const completionProvider = new completionProvider_1.AbacusCompletionProvider(dataProvider);
    const hoverProvider = new hoverProvider_1.AbacusHoverProvider(dataProvider);
    diagnosticProvider = new diagnosticProvider_1.AbacusDiagnosticProvider(dataProvider);
    const completionDisposable = vscode.languages.registerCompletionItemProvider('abacus-input', completionProvider, ...'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_ '.split(''));
    const hoverDisposable = vscode.languages.registerHoverProvider('abacus-input', hoverProvider);
    const diagnosticsDisposable = vscode.workspace.onDidChangeTextDocument((event) => {
        if (event.document.languageId === 'abacus-input' && dataLoaded) {
            const tokenSource = new vscode.CancellationTokenSource();
            diagnosticProvider?.diagnose(event.document, tokenSource.token);
        }
    });
    const openDisposable = vscode.workspace.onDidOpenTextDocument((document) => {
        if (document.languageId === 'abacus-input' && dataLoaded) {
            const tokenSource = new vscode.CancellationTokenSource();
            diagnosticProvider?.diagnose(document, tokenSource.token);
        }
    });
    const configDisposable = vscode.workspace.onDidChangeConfiguration((event) => {
        if (event.affectsConfiguration('abacusInput') && dataLoaded) {
            diagnoseAllOpenDocuments();
        }
    });
    const reloadCommand = vscode.commands.registerCommand('abacusInput.reloadData', async () => {
        try {
            await dataProvider.loadData();
            vscode.window.showInformationMessage('ABACUS INPUT data reloaded successfully');
            diagnoseAllOpenDocuments();
        }
        catch (error) {
            vscode.window.showErrorMessage(`Failed to reload data: ${error}`);
        }
    });
    const checkDataCommand = vscode.commands.registerCommand('abacusInput.checkData', async () => {
        const keywords = dataProvider.getAllKeywords();
        const completionItems = dataProvider.getCompletionItems();
        vscode.window.showInformationMessage(`ABACUS Data Status:\n` +
            `- Constraints: ${keywords.length} keywords\n` +
            `- Completion items: ${completionItems.length}\n` +
            `- Extension path: ${context.extensionPath}\n` +
            `- Data path: ${vscode.workspace.getConfiguration('abacusInput').get('dataPath', '') || '(default)'}`);
    });
    context.subscriptions.push(completionDisposable, hoverDisposable, diagnosticsDisposable, openDisposable, configDisposable, reloadCommand, checkDataCommand);
    diagnoseAllOpenDocuments();
}
function deactivate() {
    diagnosticProvider?.dispose();
    console.log('ABACUS INPUT Helper extension is now deactivated');
}
//# sourceMappingURL=extension.js.map