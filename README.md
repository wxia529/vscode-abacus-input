# ABACUS INPUT Helper for VSCode

A VSCode extension providing intelligent completion, validation, and diagnostics for ABACUS INPUT files.

## Features

- **Intelligent Completion**: Auto-completion for all ABACUS INPUT keywords with enum value dropdowns
- **Hover Documentation**: Display detailed documentation when hovering over keywords
- **Validation & Diagnostics**: 
  - Type checking for parameter values
  - Enum value validation
  - Numeric range validation
  - Dependency checking between parameters
  - Conditional warnings based on parameter combinations
- **Syntax Highlighting**: Color-coded highlighting for keywords, values, numbers, strings, and comments
- **Snippets**: Quick insertion of common parameter configurations

## Usage

1. Open any `.input` or `.INPUT` file in VSCode
2. Start typing a keyword to see completion suggestions
3. Press `Ctrl+Space` to manually trigger completion
4. Hover over any keyword to see its documentation
5. Errors and warnings will be highlighted in the editor

## Configuration

- `abacusInput.completion.enabled`: Enable/disable intelligent completion (default: true)
- `abacusInput.diagnostics.enabled`: Enable/disable validation and diagnostics (default: true)
- `abacusInput.hover.enabled`: Enable/disable hover documentation (default: true)
- `abacusInput.dataPath`: Custom path to ABACUS JSON data files (default: bundled data)

## Commands

- `ABACUS INPUT: Reload Data`: Reload the JSON data files without restarting VSCode

## Development

### Prerequisites

- Node.js >= 18
- npm >= 9
- VSCode >= 1.74

### Build

```bash
npm install
npm run compile
```

### Watch Mode

```bash
npm run watch
```

### Run Tests

```bash
npm test
```

### Package Extension

```bash
npm run package
```

### Install Extension

```bash
code --install-extension abacus-input-helper-1.0.0.vsix
```

## Project Structure

```
abacus-input/
├── src/
│   ├── extension.ts          # Main extension entry point
│   ├── dataProvider.ts       # JSON data loading and access
│   ├── completionProvider.ts # Completion provider implementation
│   ├── hoverProvider.ts      # Hover provider implementation
│   ├── diagnosticProvider.ts # Diagnostic provider implementation
│   └── test/
│       ├── extension.test.ts # Unit tests
│       └── runTest.ts        # Test runner
├── data/
│   ├── completion.json       # Completion items data
│   ├── constraints.json      # Parameter constraints
│   ├── diagnostics.json      # Diagnostic rules
│   ├── snippets.json         # Snippet definitions
│   └── dependencies.json     # Parameter dependencies
├── snippets/
│   └── abacus-input-snippets.json
├── syntaxes/
│   └── abacus-input.tmLanguage.json
├── language-configuration.json
├── package.json
└── tsconfig.json
```

## Data Files

The extension uses JSON data files generated from ABACUS documentation:

- **completion.json**: Contains all keywords with their documentation, types, and enum values
- **constraints.json**: Defines type constraints, allowed values, and ranges for each parameter
- **diagnostics.json**: Contains conditional diagnostic rules
- **snippets.json**: VSCode snippet definitions
- **dependencies.json**: Parameter dependency definitions

## License

MIT
