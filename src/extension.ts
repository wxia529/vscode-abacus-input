/**
 * ABACUS INPUT Helper Extension for VSCode
 * Provides intelligent completion, validation, and diagnostics for ABACUS INPUT files
 */

import * as vscode from 'vscode';
import { AbacusDataProvider } from './dataProvider';
import { AbacusCompletionProvider } from './completionProvider';
import { AbacusHoverProvider } from './hoverProvider';
import { AbacusDiagnosticProvider } from './diagnosticProvider';

let diagnosticProvider: AbacusDiagnosticProvider | undefined;
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

export function activate(context: vscode.ExtensionContext) {
  console.log('ABACUS INPUT Helper extension is now active');

  const dataProvider = new AbacusDataProvider(context);

  dataProvider.loadData().then(() => {
    dataLoaded = true;
    console.log('ABACUS data loaded successfully');
    diagnoseAllOpenDocuments();
  }).catch((error) => {
    console.error('Failed to load ABACUS data:', error);
    vscode.window.showErrorMessage('Failed to load ABACUS INPUT data. Please check the data path configuration.');
  });

  const completionProvider = new AbacusCompletionProvider(dataProvider);
  const hoverProvider = new AbacusHoverProvider(dataProvider);
  diagnosticProvider = new AbacusDiagnosticProvider(dataProvider);

  const completionDisposable = vscode.languages.registerCompletionItemProvider(
    'abacus-input',
    completionProvider,
    ...'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_ '.split('')
  );

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
    } catch (error) {
      vscode.window.showErrorMessage(`Failed to reload data: ${error}`);
    }
  });

  const checkDataCommand = vscode.commands.registerCommand('abacusInput.checkData', async () => {
    const keywords = dataProvider.getAllKeywords();
    const completionItems = dataProvider.getCompletionItems();
    vscode.window.showInformationMessage(
      `ABACUS Data Status:\n` +
      `- Constraints: ${keywords.length} keywords\n` +
      `- Completion items: ${completionItems.length}\n` +
      `- Extension path: ${context.extensionPath}\n` +
      `- Data path: ${vscode.workspace.getConfiguration('abacusInput').get('dataPath', '') || '(default)'}`
    );
  });

  context.subscriptions.push(
    completionDisposable,
    hoverDisposable,
    diagnosticsDisposable,
    openDisposable,
    configDisposable,
    reloadCommand,
    checkDataCommand
  );

  diagnoseAllOpenDocuments();
}

export function deactivate() {
  diagnosticProvider?.dispose();
  console.log('ABACUS INPUT Helper extension is now deactivated');
}
