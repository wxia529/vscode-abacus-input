/**
 * Data loader for ABACUS INPUT JSON files
 */

import * as vscode from 'vscode';
import * as fs from 'fs';
import * as path from 'path';

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
  secondElement?: {
    enum?: string[];
    optional?: boolean;
  };
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

export class AbacusDataProvider {
  private completionData: CompletionItem[] = [];
  private constraintsData: Record<string, Constraint> = {};
  private diagnosticsData: DiagnosticRule[] = [];
  private snippetsData: Record<string, Snippet> = {};
  private dependenciesData: Record<string, Dependency> = {};
  private dataPath: string = '';

  constructor(context: vscode.ExtensionContext) {
    this.dataPath = vscode.workspace.getConfiguration('abacusInput').get('dataPath', '') 
      || context.extensionPath;
  }

  async loadData(): Promise<void> {
    await Promise.all([
      this.loadJsonFile<CompletionItem[]>('completion.json').then(data => this.completionData = data),
      this.loadJsonFile<Record<string, Constraint>>('constraints.json').then(data => this.constraintsData = data),
      this.loadJsonFile<DiagnosticRule[]>('diagnostics.json').then(data => this.diagnosticsData = data),
      this.loadJsonFile<Record<string, Snippet>>('snippets.json').then(data => this.snippetsData = data),
      this.loadJsonFile<Record<string, Dependency>>('dependencies.json').then(data => this.dependenciesData = data)
    ]);
  }

  private async loadJsonFile<T>(filename: string): Promise<T> {
    const filePath = path.join(this.dataPath, 'data', filename);
    console.log(`Loading data from: ${filePath}`);
    try {
      const content = await fs.promises.readFile(filePath, 'utf-8');
      const data = JSON.parse(content) as T;
      console.log(`Successfully loaded ${filename}`);
      return data;
    } catch (error) {
      console.error(`Failed to load ${filename} from ${filePath}:`, error);
      vscode.window.showErrorMessage(`Failed to load ${filename}: ${error}`);
      return {} as T;
    }
  }

  getCompletionItems(): CompletionItem[] {
    return this.completionData;
  }

  getConstraint(keyword: string): Constraint | undefined {
    return this.constraintsData[keyword];
  }

  getDiagnosticRules(): DiagnosticRule[] {
    return this.diagnosticsData;
  }

  getSnippets(): Record<string, Snippet> {
    return this.snippetsData;
  }

  getDependency(keyword: string): Dependency | undefined {
    return this.dependenciesData[keyword];
  }

  getAllKeywords(): string[] {
    return Object.keys(this.constraintsData);
  }

  getCompletionByLabel(label: string): CompletionItem | undefined {
    return this.completionData.find(item => item.label === label);
  }
}
