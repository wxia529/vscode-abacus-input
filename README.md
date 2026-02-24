# ABACUS INPUT for VS Code

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Keyword completion, validation, and diagnostics for ABACUS INPUT files in Visual Studio Code.

## Features

This extension provides support for editing ABACUS INPUT files:

- **🎯 Keyword Completion** - Keyword and value completion based on ABACUS documentation
- **✅ Validation & Diagnostics** - Error detection and parameter validation
- **📖 Hover Documentation** - Keyword documentation on hover
- **🎨 Syntax Highlighting** - Syntax highlighting for INPUT files
- **📝 Code Snippets** - Quick insertion for common ABACUS parameters

## Installation

Install from the [VS Code Marketplace](https://marketplace.visualstudio.com/items?itemName=wxia529.abacus-input) or search for "ABACUS INPUT" in the VS Code Extensions panel.

## Usage

### Automatic Activation

The extension automatically activates when:
- Opening a file named `INPUT` or `INPUT.*`
- A file with ABACUS INPUT content is opened

### Commands

| Command | Description |
|---------|-------------|
| `ABACUS INPUT: Reload Data` | Reload the completion data from JSON files |
| `ABACUS INPUT: Check Data Status` | Display current data status and statistics |

### Configuration

| Setting | Default | Description |
|---------|---------|-------------|
| `abacusInput.completion.enabled` | `true` | Enable completion |
| `abacusInput.diagnostics.enabled` | `true` | Enable validation and diagnostics |
| `abacusInput.hover.enabled` | `true` | Enable hover documentation |
| `abacusInput.dataPath` | `""` | Custom path to ABACUS JSON data files |

## Development

### Build

```bash
npm install
npm run compile
```

### Watch Mode

```bash
npm run watch
```

### Lint

```bash
npm run lint
```

### Package

```bash
npm run package
```

## Project Structure

```
vscode-abacus-input/
├── src/                    # Extension source code
│   ├── extension.ts        # Main extension entry point
│   ├── completionProvider.ts
│   ├── diagnosticProvider.ts
│   ├── hoverProvider.ts
│   └── dataProvider.ts
├── syntaxes/               # TextMate grammar for syntax highlighting
├── snippets/               # Code snippets
├── data/                   # JSON data files for completion
├── language-configuration.json
└── package.json            # Extension manifest
```

## License

MIT License - see [LICENSE](LICENSE) for details.

## Acknowledgments

- [ABACUS](https://github.com/deepmodeling/abacus-develop) - Atomic-orbital Based Ab-initio Computation at UStc
- Data sourced from the official ABACUS documentation
