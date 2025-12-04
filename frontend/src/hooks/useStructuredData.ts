/**
 * =============================================================================
 * USE STRUCTURED DATA HOOK
 * =============================================================================
 * 
 * Hook for MainPage to subscribe to structured data events from chat responses.
 * 
 * When the chat receives a response with structured data (molecule, bioisosteres, etc.),
 * it dispatches a custom event. This hook lets components subscribe to those events.
 * 
 * USAGE:
 * ```typescript
 * const { latestData, clearData } = useStructuredData();
 * 
 * useEffect(() => {
 *   if (latestData?.action === 'RENDER_MOLECULE' && latestData.moleculeData) {
 *     setMolecule(latestData.moleculeData);
 *   }
 * }, [latestData]);
 * ```
 */

import { useState, useEffect, useCallback } from 'react';
import type { ParsedChatResponse, RenderAction } from '@/utils/chatResponseHandler';

// =============================================================================
// EVENT NAME CONSTANTS
// =============================================================================

/** Custom event name for structured data dispatch */
export const STRUCTURED_DATA_EVENT = 'sedona:structured-data';

/** Custom event name for clearing data */
export const CLEAR_STRUCTURED_DATA_EVENT = 'sedona:clear-structured-data';

// =============================================================================
// EVENT TYPES
// =============================================================================

export interface StructuredDataEvent extends CustomEvent {
  detail: ParsedChatResponse;
}

// =============================================================================
// EVENT DISPATCH FUNCTIONS (used by ChatContext)
// =============================================================================

/**
 * Dispatch a structured data event.
 * Called by ChatContext when a response contains structured data.
 * 
 * @param data - The parsed chat response with structured data
 */
export function dispatchStructuredData(data: ParsedChatResponse): void {
  const event = new CustomEvent(STRUCTURED_DATA_EVENT, {
    detail: data,
    bubbles: true,
  });
  window.dispatchEvent(event);
  console.log('📡 Dispatched structured data event:', data.action);
}

/**
 * Dispatch an event to clear structured data.
 * Can be used to reset the MainPage visualization.
 */
export function dispatchClearStructuredData(): void {
  const event = new CustomEvent(CLEAR_STRUCTURED_DATA_EVENT);
  window.dispatchEvent(event);
  console.log('📡 Dispatched clear structured data event');
}

// =============================================================================
// HOOK OPTIONS
// =============================================================================

export interface UseStructuredDataOptions {
  /** Only subscribe to specific action types */
  filterActions?: RenderAction[];
  
  /** Auto-clear data after a delay (ms) - 0 means never auto-clear */
  autoClearDelay?: number;
  
  /** Called when new data arrives */
  onData?: (data: ParsedChatResponse) => void;
}

// =============================================================================
// MAIN HOOK
// =============================================================================

/**
 * Subscribe to structured data events from chat responses.
 * 
 * @param options - Configuration options
 * @returns Object with current data and control functions
 * 
 * @example
 * ```typescript
 * // Subscribe to all structured data
 * const { latestData } = useStructuredData();
 * 
 * // Subscribe only to molecule and bioisostere data
 * const { latestData } = useStructuredData({
 *   filterActions: ['RENDER_MOLECULE', 'RENDER_BIOISOSTERES'],
 * });
 * 
 * // With callback
 * useStructuredData({
 *   onData: (data) => {
 *     if (data.action === 'RENDER_MOLECULE') {
 *       handleMolecule(data.moleculeData);
 *     }
 *   }
 * });
 * ```
 */
export function useStructuredData(options: UseStructuredDataOptions = {}) {
  const { filterActions, autoClearDelay = 0, onData } = options;
  
  const [latestData, setLatestData] = useState<ParsedChatResponse | null>(null);
  const [history, setHistory] = useState<ParsedChatResponse[]>([]);
  
  // Clear the latest data
  const clearData = useCallback(() => {
    setLatestData(null);
  }, []);
  
  // Clear all history
  const clearHistory = useCallback(() => {
    setHistory([]);
    setLatestData(null);
  }, []);
  
  // Handle incoming structured data
  useEffect(() => {
    const handleStructuredData = (event: Event) => {
      const customEvent = event as StructuredDataEvent;
      const data = customEvent.detail;
      
      // Apply filter if specified
      if (filterActions && filterActions.length > 0) {
        if (!filterActions.includes(data.action)) {
          return; // Skip this event
        }
      }
      
      // Update state
      setLatestData(data);
      setHistory(prev => [...prev, data]);
      
      // Call callback if provided
      if (onData) {
        onData(data);
      }
      
      // Auto-clear after delay if specified
      if (autoClearDelay > 0) {
        setTimeout(() => {
          setLatestData(null);
        }, autoClearDelay);
      }
    };
    
    const handleClearData = () => {
      clearData();
    };
    
    // Subscribe to events
    window.addEventListener(STRUCTURED_DATA_EVENT, handleStructuredData);
    window.addEventListener(CLEAR_STRUCTURED_DATA_EVENT, handleClearData);
    
    // Cleanup
    return () => {
      window.removeEventListener(STRUCTURED_DATA_EVENT, handleStructuredData);
      window.removeEventListener(CLEAR_STRUCTURED_DATA_EVENT, handleClearData);
    };
  }, [filterActions, autoClearDelay, onData, clearData]);
  
  return {
    /** The most recent structured data received */
    latestData,
    
    /** History of all structured data received */
    history,
    
    /** Clear the latest data */
    clearData,
    
    /** Clear all history */
    clearHistory,
    
    /** Whether there is currently data to display */
    hasData: latestData !== null,
  };
}

// =============================================================================
// SPECIALIZED HOOKS (convenience wrappers)
// =============================================================================

/**
 * Subscribe only to molecule annotation data.
 */
export function useMoleculeData() {
  return useStructuredData({
    filterActions: ['RENDER_MOLECULE'],
  });
}

/**
 * Subscribe only to bioisostere data.
 */
export function useBioisostereData() {
  return useStructuredData({
    filterActions: ['RENDER_BIOISOSTERES', 'RENDER_MULTI_BIO'],
  });
}

/**
 * Subscribe only to ADMET data.
 */
export function useAdmetData() {
  return useStructuredData({
    filterActions: ['RENDER_ADMET', 'RENDER_MULTI_ADMET'],
  });
}

/**
 * Subscribe only to query set data.
 */
export function useQuerySetData() {
  return useStructuredData({
    filterActions: ['RENDER_QUERY_SET'],
  });
}

// =============================================================================
// MOLECULE CONTEXT STORAGE (for conversation history)
// =============================================================================

/**
 * Pattern info for a detected substructure.
 */
export interface PatternInfo {
  letter: string;       // A, B, C, etc.
  name: string;         // "Pyridine", "Piperazine", etc.
  smiles: string;       // SMILES/SMARTS of the pattern
  atomIndices: number[]; // Atom indices in parent molecule
}

/**
 * Current molecule context that can be passed to the backend.
 */
export interface MoleculeContext {
  /** The full molecule SMILES */
  smiles: string;
  /** Common name of the molecule (e.g., "Sildenafil") */
  name?: string;
  /** Detected patterns mapped by letter (A, B, C...) */
  patterns: PatternInfo[];
  /** Timestamp when context was set */
  timestamp: number;
}

// Module-level storage for current molecule context
let _currentMoleculeContext: MoleculeContext | null = null;

/**
 * Set the current molecule context.
 * Called by MainPage when molecule analysis changes.
 * 
 * @param context - The molecule context or null to clear
 */
export function setMoleculeContext(context: MoleculeContext | null): void {
  _currentMoleculeContext = context;
  if (context) {
    console.log('🧬 Molecule context updated:', {
      smiles: context.smiles.substring(0, 30) + '...',
      patterns: context.patterns.length
    });
  } else {
    console.log('🧬 Molecule context cleared');
  }
}

/**
 * Get the current molecule context.
 * Called by ChatContext when building conversation history.
 * 
 * @returns The current molecule context or null
 */
export function getMoleculeContext(): MoleculeContext | null {
  return _currentMoleculeContext;
}

/**
 * Get a formatted JSON string of the molecule context for conversation history.
 * This is what gets appended to messages for the LLM to understand.
 * 
 * @returns A formatted string describing the current molecule context, or empty string if none
 */
export function getMoleculeContextString(): string {
  if (!_currentMoleculeContext) return '';
  
  const ctx = _currentMoleculeContext;
  
  // Build a compact but informative context string
  const patternDescriptions = ctx.patterns.map(p => 
    `  ${p.letter}: "${p.name}" (${p.smiles})`
  ).join('\n');
  
  const moleculeName = ctx.name && ctx.name !== 'No Name' && ctx.name !== 'Unknown'
    ? `Molecule Name: ${ctx.name}\n`
    : '';
  
  return `[CURRENT MOLECULE CONTEXT]
${moleculeName}Molecule SMILES: ${ctx.smiles}
Detected Substructures:
${patternDescriptions || '  (none detected)'}
[END CONTEXT]`;
}

/**
 * Clear the molecule context.
 */
export function clearMoleculeContext(): void {
  _currentMoleculeContext = null;
  console.log('🧬 Molecule context cleared');
}
