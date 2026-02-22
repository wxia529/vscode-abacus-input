"use strict";
/**
 * Data loader for ABACUS INPUT JSON files
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
exports.AbacusDataProvider = void 0;
const vscode = __importStar(require("vscode"));
const fs = __importStar(require("fs"));
const path = __importStar(require("path"));
class AbacusDataProvider {
    constructor(context) {
        this.completionData = [];
        this.constraintsData = {};
        this.diagnosticsData = [];
        this.snippetsData = {};
        this.dependenciesData = {};
        this.dataPath = '';
        this.dataPath = vscode.workspace.getConfiguration('abacusInput').get('dataPath', '')
            || context.extensionPath;
    }
    async loadData() {
        await Promise.all([
            this.loadJsonFile('completion.json').then(data => this.completionData = data),
            this.loadJsonFile('constraints.json').then(data => this.constraintsData = data),
            this.loadJsonFile('diagnostics.json').then(data => this.diagnosticsData = data),
            this.loadJsonFile('snippets.json').then(data => this.snippetsData = data),
            this.loadJsonFile('dependencies.json').then(data => this.dependenciesData = data)
        ]);
    }
    async loadJsonFile(filename) {
        const filePath = path.join(this.dataPath, 'data', filename);
        console.log(`Loading data from: ${filePath}`);
        try {
            const content = await fs.promises.readFile(filePath, 'utf-8');
            const data = JSON.parse(content);
            console.log(`Successfully loaded ${filename}`);
            return data;
        }
        catch (error) {
            console.error(`Failed to load ${filename} from ${filePath}:`, error);
            vscode.window.showErrorMessage(`Failed to load ${filename}: ${error}`);
            return {};
        }
    }
    getCompletionItems() {
        return this.completionData;
    }
    getConstraint(keyword) {
        return this.constraintsData[keyword];
    }
    getDiagnosticRules() {
        return this.diagnosticsData;
    }
    getSnippets() {
        return this.snippetsData;
    }
    getDependency(keyword) {
        return this.dependenciesData[keyword];
    }
    getAllKeywords() {
        return Object.keys(this.constraintsData);
    }
    getCompletionByLabel(label) {
        return this.completionData.find(item => item.label === label);
    }
}
exports.AbacusDataProvider = AbacusDataProvider;
//# sourceMappingURL=dataProvider.js.map