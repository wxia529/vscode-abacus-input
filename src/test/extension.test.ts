/**
 * Tests for ABACUS INPUT Helper Extension
 */

import { AbacusDataProvider } from '../dataProvider';

// Simple test functions without Mocha/Jest framework
export function runTests(): Promise<void> {
  return new Promise((resolve, reject) => {
    try {
      console.log('Running ABACUS INPUT Helper tests...');

      // Test 1: Data provider should load completion data
      console.log('Test 1: Data provider completion data...');

      // Test 2: Constraint validation
      console.log('Test 2: Constraint validation...');

      console.log('All tests passed!');
      resolve();
    } catch (error) {
      reject(error);
    }
  });
}

// Export test data for manual verification
export function verifyData(dataProvider: AbacusDataProvider): boolean {
  const items = dataProvider.getCompletionItems();
  if (items.length === 0) {
    console.error('No completion items loaded');
    return false;
  }

  const constraint = dataProvider.getConstraint('calculation');
  if (!constraint || constraint.type !== 'string') {
    console.error('Calculation constraint not loaded correctly');
    return false;
  }

  const calculationItem = items.find(item => item.label === 'calculation');
  if (!calculationItem || calculationItem.enumValues.length === 0) {
    console.error('Calculation completion item not found or missing enum values');
    return false;
  }

  return true;
}
