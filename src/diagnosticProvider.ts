/**
 * Diagnostic Provider for ABACUS INPUT files
 * Validates parameters against constraints and dependency rules
 */

import * as vscode from 'vscode';
import { AbacusDataProvider, Constraint } from './dataProvider';

export class AbacusDiagnosticProvider {
  private dataProvider: AbacusDataProvider;
  private diagnosticCollection: vscode.DiagnosticCollection;

  constructor(dataProvider: AbacusDataProvider) {
    this.dataProvider = dataProvider;
    this.diagnosticCollection = vscode.languages.createDiagnosticCollection('abacus-input');
  }

  public dispose(): void {
    this.diagnosticCollection.dispose();
  }

  public async diagnose(
    document: vscode.TextDocument,
    token: vscode.CancellationToken
  ): Promise<void> {
    const config = vscode.workspace.getConfiguration('abacusInput');
    if (!config.get('diagnostics.enabled', true)) {
      this.diagnosticCollection.delete(document.uri);
      return;
    }

    if (document.languageId !== 'abacus-input') {
      return;
    }

    const diagnostics: vscode.Diagnostic[] = [];
    const { keywordValues, duplicateDiagnostics } = this.parseDocument(document);
    diagnostics.push(...duplicateDiagnostics);

    for (let lineNum = 0; lineNum < document.lineCount; lineNum++) {
      if (token.isCancellationRequested) return;
      const line = document.lineAt(lineNum);
      diagnostics.push(...this.validateLine(line, keywordValues));
    }

    this.diagnosticCollection.set(document.uri, diagnostics);
  }

  private parseDocument(document: vscode.TextDocument): {
    keywordValues: Map<string, { value: string; line: number }>;
    duplicateDiagnostics: vscode.Diagnostic[];
  } {
    const keywordValues = new Map<string, { value: string; line: number }>();
    const duplicateKeywords = new Map<string, number[]>();
    const duplicateDiagnostics: vscode.Diagnostic[] = [];

    for (let i = 0; i < document.lineCount; i++) {
      const line = document.lineAt(i);
      const match = this.parseLine(line.text);
      if (match) {
        const { keyword } = match;

        if (duplicateKeywords.has(keyword)) {
          duplicateKeywords.get(keyword)!.push(i);
        } else if (keywordValues.has(keyword)) {
          const originalLine = keywordValues.get(keyword)!.line;
          duplicateKeywords.set(keyword, [originalLine, i]);
          keywordValues.delete(keyword);
        } else {
          keywordValues.set(keyword, { value: match.value, line: i });
        }
      }
    }

    for (const [keyword, lines] of duplicateKeywords.entries()) {
      for (const lineNum of lines) {
        const line = document.lineAt(lineNum);
        const diagnostic = new vscode.Diagnostic(
          line.range,
          `Duplicate keyword: '${keyword}' is defined multiple times`,
          vscode.DiagnosticSeverity.Warning
        );
        diagnostic.code = 'duplicate-keyword';
        diagnostic.source = 'abacus-input';
        duplicateDiagnostics.push(diagnostic);
      }
    }

    return { keywordValues, duplicateDiagnostics };
  }

  private parseLine(text: string): { keyword: string; value: string } | null {
    const trimmed = text.trim();
    if (!trimmed || trimmed.startsWith('#')) return null;

    const match = trimmed.match(/^(\w+)(?:\s+(.*))?$/);
    if (match) {
      let value = match[2] ? match[2].trim() : '';
      const commentIndex = value.indexOf('#');
      if (commentIndex !== -1) value = value.substring(0, commentIndex).trim();
      return { keyword: match[1], value };
    }

    return null;
  }

  private validateLine(
    line: vscode.TextLine,
    keywordValues: Map<string, { value: string; line: number }>
  ): vscode.Diagnostic[] {
    const diagnostics: vscode.Diagnostic[] = [];
    const match = this.parseLine(line.text);
    if (!match) return diagnostics;

    const { keyword, value } = match;
    const constraint = this.dataProvider.getConstraint(keyword);

    // Check for unknown keywords
    if (!constraint) {
      const allowedKeywords = ['INPUT_PARAMETERS'];
      if (!line.text.trim().startsWith('#') && !allowedKeywords.includes(keyword)) {
        diagnostics.push(this.createDiagnostic(
          line.range,
          `Unknown keyword: '${keyword}'`,
          vscode.DiagnosticSeverity.Warning,
          'unknown-keyword'
        ));
      }
      return diagnostics;
    }

    // Check for missing value
    if (!value) {
      diagnostics.push(this.createDiagnostic(
        line.range,
        `Missing value for keyword: '${keyword}'`,
        vscode.DiagnosticSeverity.Error,
        'missing-value'
      ));
      return diagnostics;
    }

    diagnostics.push(
      ...this.validateType(keyword, value, constraint, line.range),
      ...this.validateEnum(keyword, value, constraint, line.range),
      ...this.validateRange(keyword, value, constraint, line.range),
      ...this.validateDependencies(keyword, value, keywordValues, line.range),
      ...this.validateConditionalRules(keyword, value, keywordValues, line.range)
    );

    return diagnostics;
  }

  private validateType(
    keyword: string,
    value: string,
    constraint: Constraint,
    range: vscode.Range
  ): vscode.Diagnostic[] {
    const expectedType = constraint.type;
    const diagnostics: vscode.Diagnostic[] = [];

    if (expectedType === 'string') return diagnostics;

    if (expectedType === 'number') {
      if (isNaN(parseFloat(value))) {
        diagnostics.push(this.createDiagnostic(
          range,
          `Expected a number for '${keyword}', got '${value}'`,
          vscode.DiagnosticSeverity.Error,
          'type-mismatch'
        ));
      }
      return diagnostics;
    }

    if (expectedType === 'integer') {
      if (isNaN(parseInt(value, 10)) || !/^-?\d+$/.test(value.trim())) {
        diagnostics.push(this.createDiagnostic(
          range,
          `Expected an integer for '${keyword}', got '${value}'`,
          vscode.DiagnosticSeverity.Error,
          'type-mismatch'
        ));
      }
      return diagnostics;
    }

    if (expectedType === 'boolean') {
      const validBooleans = ['true', 'false', 'TRUE', 'FALSE', 'True', 'False', '1', '0'];
      if (!validBooleans.includes(value.trim())) {
        diagnostics.push(this.createDiagnostic(
          range,
          `Expected a boolean for '${keyword}', got '${value}'`,
          vscode.DiagnosticSeverity.Error,
          'type-mismatch'
        ));
      }
    }

    return diagnostics;
  }

  private validateEnum(
    keyword: string,
    value: string,
    constraint: Constraint,
    range: vscode.Range
  ): vscode.Diagnostic[] {
    if (!constraint.enum || constraint.enum.length === 0) return [];

    if (!constraint.enum.includes(value)) {
      return [this.createDiagnostic(
        range,
        `Invalid value '${value}' for '${keyword}'. Allowed values: ${constraint.enum.join(', ')}`,
        vscode.DiagnosticSeverity.Error,
        'invalid-enum'
      )];
    }

    return [];
  }

  private validateRange(
    keyword: string,
    value: string,
    constraint: Constraint,
    range: vscode.Range
  ): vscode.Diagnostic[] {
    const numValue = parseFloat(value);
    if (isNaN(numValue)) return [];

    const diagnostics: vscode.Diagnostic[] = [];

    if (constraint.minimum !== undefined && numValue < constraint.minimum) {
      diagnostics.push(this.createDiagnostic(
        range,
        `Value ${numValue} is less than minimum ${constraint.minimum} for '${keyword}'`,
        vscode.DiagnosticSeverity.Error,
        'range-violation'
      ));
    }

    if (constraint.maximum !== undefined && numValue > constraint.maximum) {
      diagnostics.push(this.createDiagnostic(
        range,
        `Value ${numValue} is greater than maximum ${constraint.maximum} for '${keyword}'`,
        vscode.DiagnosticSeverity.Error,
        'range-violation'
      ));
    }

    if (constraint.exclusiveMinimum !== undefined && numValue <= constraint.exclusiveMinimum) {
      diagnostics.push(this.createDiagnostic(
        range,
        `Value ${numValue} must be greater than ${constraint.exclusiveMinimum} for '${keyword}'`,
        vscode.DiagnosticSeverity.Error,
        'range-violation'
      ));
    }

    if (constraint.exclusiveMaximum !== undefined && numValue >= constraint.exclusiveMaximum) {
      diagnostics.push(this.createDiagnostic(
        range,
        `Value ${numValue} must be less than ${constraint.exclusiveMaximum} for '${keyword}'`,
        vscode.DiagnosticSeverity.Error,
        'range-violation'
      ));
    }

    return diagnostics;
  }

  private validateDependencies(
    keyword: string,
    value: string,
    keywordValues: Map<string, { value: string; line: number }>,
    range: vscode.Range
  ): vscode.Diagnostic[] {
    const diagnostics: vscode.Diagnostic[] = [];
    const dependency = this.dataProvider.getDependency(keyword);
    if (!dependency) return diagnostics;

    for (const [requiredKeyword, requiredValue] of Object.entries(dependency.requires)) {
      const actualValue = keywordValues.get(requiredKeyword);
      const constraint = this.dataProvider.getConstraint(requiredKeyword);

      if (!actualValue) {
        if (constraint?.default !== undefined && this.valuesMatch(constraint.default, requiredValue)) {
          continue;
        }
        diagnostics.push(this.createDiagnostic(
          range,
          `'${keyword}' requires '${requiredKeyword}' to be set`,
          vscode.DiagnosticSeverity.Warning,
          'missing-dependency'
        ));
        continue;
      }

      if (!this.valuesMatch(actualValue.value, requiredValue)) {
        diagnostics.push(this.createDiagnostic(
          range,
          `'${keyword}' requires '${requiredKeyword}' to be '${requiredValue}', but it is '${actualValue.value}'`,
          vscode.DiagnosticSeverity.Warning,
          'dependency-mismatch'
        ));
      }
    }

    return diagnostics;
  }

  private validateConditionalRules(
    keyword: string,
    value: string,
    keywordValues: Map<string, { value: string; line: number }>,
    range: vscode.Range
  ): vscode.Diagnostic[] {
    const diagnostics: vscode.Diagnostic[] = [];
    const rules = this.dataProvider.getDiagnosticRules();

    for (const rule of rules) {
      if (rule.keyword !== keyword) continue;

      if (rule.condition) {
        const { keyword: condKeyword, operator, value: condValue } = rule.condition;
        const actualValue = keywordValues.get(condKeyword);
        if (!actualValue || !this.evaluateCondition(actualValue.value, operator, condValue)) {
          continue;
        }
      }

      const diagnostic = new vscode.Diagnostic(range, rule.message, this.getDiagnosticSeverity(rule.severity));
      diagnostic.code = 'conditional-rule';
      diagnostic.source = 'abacus-input';
      diagnostics.push(diagnostic);
    }

    return diagnostics;
  }

  private evaluateCondition(
    actualValue: string,
    operator: string,
    expectedValue: string | number | boolean | string[]
  ): boolean {
    const actual = actualValue.trim();

    switch (operator) {
      case '==': return actual === String(expectedValue);
      case '!=': return actual !== String(expectedValue);
      case 'in': return Array.isArray(expectedValue) && expectedValue.includes(actual);
      case '>': return parseFloat(actual) > Number(expectedValue);
      case '>=': return parseFloat(actual) >= Number(expectedValue);
      case '<': return parseFloat(actual) < Number(expectedValue);
      case '<=': return parseFloat(actual) <= Number(expectedValue);
      default: return false;
    }
  }

  private valuesMatch(actual: string, expected: string | number | boolean): boolean {
    const actualTrimmed = actual.trim();
    const expectedStr = String(expected);

    if (actualTrimmed === expectedStr) return true;

    // Boolean normalization
    const boolMap: Record<string, boolean> = {
      'true': true, 'TRUE': true, 'True': true, '1': true,
      'false': false, 'FALSE': false, 'False': false, '0': false
    };

    if (typeof expected === 'boolean') {
      return boolMap[actualTrimmed] === expected;
    }

    // Number comparison
    const actualNum = parseFloat(actualTrimmed);
    const expectedNum = Number(expected);
    if (!isNaN(actualNum) && !isNaN(expectedNum)) {
      return actualNum === expectedNum;
    }

    return false;
  }

  private getDiagnosticSeverity(severity: string): vscode.DiagnosticSeverity {
    switch (severity.toLowerCase()) {
      case 'error': return vscode.DiagnosticSeverity.Error;
      case 'warning': return vscode.DiagnosticSeverity.Warning;
      default: return vscode.DiagnosticSeverity.Information;
    }
  }

  private createDiagnostic(
    range: vscode.Range,
    message: string,
    severity: vscode.DiagnosticSeverity,
    code: string
  ): vscode.Diagnostic {
    const diagnostic = new vscode.Diagnostic(range, message, severity);
    diagnostic.code = code;
    diagnostic.source = 'abacus-input';
    return diagnostic;
  }
}
