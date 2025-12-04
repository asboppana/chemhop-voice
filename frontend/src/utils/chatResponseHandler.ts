/**
 * =============================================================================
 * CHAT RESPONSE HANDLER
 * =============================================================================
 * 
 * This module maps backend AgentChatResponse types to frontend rendering actions.
 * 
 * The backend (voice_agent.py) sends responses with these possible types:
 * 
 * ┌─────────────────────┬─────────────────────────────────────────────────────┐
 * │ response_type       │ Description                                         │
 * ├─────────────────────┼─────────────────────────────────────────────────────┤
 * │ text_only           │ Pure text response - render in chat bubble          │
 * │ annotation          │ Molecule with SVG + functional groups               │
 * │ query_set           │ List of SMILES targets to query                     │
 * │ bioisosteres        │ Single bioisostere scan result                      │
 * │ admet               │ Single ADMET property predictions                   │
 * │ multi_bioisostere   │ Multiple bioisostere results                        │
 * │ multi_admet         │ Multiple ADMET predictions                          │
 * └─────────────────────┴─────────────────────────────────────────────────────┘
 * 
 * USAGE:
 * ```typescript
 * import { parseAgentResponse, isStructuredResponse, type ParsedChatResponse } from '@/utils/chatResponseHandler';
 * 
 * const response = await chatWithAgent(messages);
 * const parsed = parseAgentResponse(response);
 * 
 * if (parsed.hasStructuredData) {
 *   // Dispatch to MainPage for visualization
 *   dispatchStructuredData(parsed);
 * }
 * ```
 */

import type {
  AgentChatResponse,
  MoleculeAnalysisResponse,
  SmilesQuerySet,
  BioisostereScanResult,
  ADMETPropertyResult,
  ChatResponseType,
} from '@/services/moleculeService';

// =============================================================================
// TYPE DEFINITIONS
// =============================================================================

/**
 * Rendering action types for the frontend.
 * Each action corresponds to a specific visualization component.
 */
export type RenderAction =
  | 'RENDER_TEXT'           // Show text in chat bubble only
  | 'RENDER_MOLECULE'       // Show molecule SVG + functional groups in MainPage
  | 'RENDER_QUERY_SET'      // Show query targets (before executing queries)
  | 'RENDER_BIOISOSTERES'   // Show bioisostere scan results
  | 'RENDER_ADMET'          // Show ADMET property predictions
  | 'RENDER_MULTI_BIO'      // Show multiple bioisostere results
  | 'RENDER_MULTI_ADMET';   // Show multiple ADMET results

/**
 * Parsed response with typed data ready for rendering.
 * Only ONE of the data fields will be populated based on responseType.
 */
export interface ParsedChatResponse {
  /** Whether this response contains structured data for visualization */
  hasStructuredData: boolean;
  
  /** The rendering action to take */
  action: RenderAction;
  
  /** Original response type from backend */
  responseType: ChatResponseType;
  
  /** Text content for chat display (always present) */
  textContent: string;
  
  // ─────────────────────────────────────────────────────────────────────────
  // TYPED DATA FIELDS (mutually exclusive - only one will be set)
  // ─────────────────────────────────────────────────────────────────────────
  
  /** Molecule annotation data (when action === 'RENDER_MOLECULE') */
  moleculeData?: MoleculeAnalysisResponse;
  
  /** Query set data (when action === 'RENDER_QUERY_SET') */
  querySetData?: SmilesQuerySet;
  
  /** Bioisostere scan result (when action === 'RENDER_BIOISOSTERES') */
  bioisostereData?: BioisostereScanResult;
  
  /** ADMET predictions (when action === 'RENDER_ADMET') */
  admetData?: ADMETPropertyResult;
  
  /** Multiple bioisostere results (when action === 'RENDER_MULTI_BIO') */
  multiBioisostereData?: BioisostereScanResult[];
  
  /** Multiple ADMET results (when action === 'RENDER_MULTI_ADMET') */
  multiAdmetData?: ADMETPropertyResult[];
}

// =============================================================================
// RESPONSE TYPE TO ACTION MAPPING
// =============================================================================

/**
 * Maps backend response_type to frontend RenderAction.
 * 
 * This is the core mapping table that determines what visualization
 * to render based on the response type.
 */
const RESPONSE_TYPE_TO_ACTION: Record<ChatResponseType, RenderAction> = {
  'text_only':         'RENDER_TEXT',
  'annotation':        'RENDER_MOLECULE',
  'query_set':         'RENDER_QUERY_SET',
  'bioisosteres':      'RENDER_BIOISOSTERES',
  'admet':             'RENDER_ADMET',
  'multi_bioisostere': 'RENDER_MULTI_BIO',
  'multi_admet':       'RENDER_MULTI_ADMET',
};

/**
 * Response types that contain structured data requiring visualization
 * beyond simple text rendering.
 */
const STRUCTURED_RESPONSE_TYPES: ChatResponseType[] = [
  'annotation',
  'query_set',
  'bioisosteres',
  'admet',
  'multi_bioisostere',
  'multi_admet',
];

// =============================================================================
// MAIN PARSING FUNCTION
// =============================================================================

/**
 * Parse an AgentChatResponse into a typed ParsedChatResponse.
 * 
 * This function:
 * 1. Determines the response type and corresponding render action
 * 2. Extracts typed data from the appropriate field
 * 3. Returns a clean interface for rendering components
 * 
 * @param response - Raw response from the backend chat API
 * @returns ParsedChatResponse with typed data ready for rendering
 * 
 * @example
 * ```typescript
 * const response = await chatWithAgent(messages);
 * const parsed = parseAgentResponse(response);
 * 
 * switch (parsed.action) {
 *   case 'RENDER_MOLECULE':
 *     renderMolecule(parsed.moleculeData!);
 *     break;
 *   case 'RENDER_BIOISOSTERES':
 *     renderBioisosteres(parsed.bioisostereData!);
 *     break;
 *   // ...
 * }
 * ```
 */
export function parseAgentResponse(response: AgentChatResponse): ParsedChatResponse {
  const responseType: ChatResponseType = response.response_type || 'text_only';
  const action = RESPONSE_TYPE_TO_ACTION[responseType] || 'RENDER_TEXT';
  const hasStructuredData = STRUCTURED_RESPONSE_TYPES.includes(responseType);
  
  const parsed: ParsedChatResponse = {
    hasStructuredData,
    action,
    responseType,
    textContent: response.response || '',
  };
  
  // Populate typed data based on response type
  switch (responseType) {
    case 'annotation':
      parsed.moleculeData = response.molecule_data || extractMoleculeFromStructured(response);
      break;
      
    case 'query_set':
      parsed.querySetData = response.query_set_data || extractQuerySetFromStructured(response);
      break;
      
    case 'bioisosteres':
      parsed.bioisostereData = response.bioisostere_data || extractBioisostereFromStructured(response);
      break;
      
    case 'admet':
      parsed.admetData = response.admet_data || extractAdmetFromStructured(response);
      break;
      
    case 'multi_bioisostere':
      parsed.multiBioisostereData = extractMultiBioisostereFromStructured(response);
      break;
      
    case 'multi_admet':
      parsed.multiAdmetData = extractMultiAdmetFromStructured(response);
      break;
  }
  
  return parsed;
}

// =============================================================================
// HELPER FUNCTIONS FOR EXTRACTING DATA FROM STRUCTURED_DATA
// =============================================================================

/**
 * Extract molecule annotation from structured_data.embedded_objects.
 * Used as fallback when molecule_data is not populated.
 */
function extractMoleculeFromStructured(response: AgentChatResponse): MoleculeAnalysisResponse | undefined {
  const objects = response.structured_data?.embedded_objects;
  if (!objects?.length) return undefined;
  
  const moleculeObj = objects.find(obj => obj.type === 'molecule_structured');
  if (!moleculeObj) return undefined;
  
  // Handle both formats: { annotation: {...} } and direct fields
  if (moleculeObj.annotation) {
    return moleculeObj.annotation as MoleculeAnalysisResponse;
  }
  
  // Direct format (less common)
  if (moleculeObj.smiles && moleculeObj.svg) {
    return {
      name: moleculeObj.name || 'Unknown',
      smiles: moleculeObj.smiles,
      svg: moleculeObj.svg,
      matches: moleculeObj.matches || [],
    };
  }
  
  return undefined;
}

/**
 * Extract query set from structured_data.embedded_objects.
 */
function extractQuerySetFromStructured(response: AgentChatResponse): SmilesQuerySet | undefined {
  const objects = response.structured_data?.embedded_objects;
  if (!objects?.length) return undefined;
  
  const querySetObj = objects.find(obj => obj.type === 'query_set_structured');
  if (!querySetObj) return undefined;
  
  if (querySetObj.query_set) {
    return querySetObj.query_set as SmilesQuerySet;
  }
  
  // Handle direct format
  if (querySetObj.targets) {
    return querySetObj as SmilesQuerySet;
  }
  
  return undefined;
}

/**
 * Extract single bioisostere result from structured_data.embedded_objects.
 */
function extractBioisostereFromStructured(response: AgentChatResponse): BioisostereScanResult | undefined {
  const objects = response.structured_data?.embedded_objects;
  if (!objects?.length) return undefined;
  
  const bioObj = objects.find(obj => obj.type === 'bioisostere_structured');
  if (!bioObj) return undefined;
  
  return {
    query_smiles: bioObj.query_smiles || '',
    source_group: bioObj.source_group,
    atom_indices: bioObj.atom_indices,
    num_results: bioObj.results?.length || 0,
    results: bioObj.results || [],
  };
}

/**
 * Extract ADMET predictions from structured_data.embedded_objects.
 */
function extractAdmetFromStructured(response: AgentChatResponse): ADMETPropertyResult | undefined {
  const objects = response.structured_data?.embedded_objects;
  if (!objects?.length) return undefined;
  
  const admetObj = objects.find(obj => obj.type === 'admet_structured');
  if (!admetObj) return undefined;
  
  return {
    smiles: admetObj.smiles || '',
    predictions: admetObj.predictions || {},
  };
}

/**
 * Extract multiple bioisostere results from structured_data.embedded_objects.
 */
function extractMultiBioisostereFromStructured(response: AgentChatResponse): BioisostereScanResult[] | undefined {
  const objects = response.structured_data?.embedded_objects;
  if (!objects?.length) return undefined;
  
  const bioObjects = objects.filter(obj => obj.type === 'bioisostere_structured');
  if (!bioObjects.length) return undefined;
  
  return bioObjects.map(bioObj => ({
    query_smiles: bioObj.query_smiles || '',
    source_group: bioObj.source_group,
    atom_indices: bioObj.atom_indices,
    num_results: bioObj.results?.length || 0,
    results: bioObj.results || [],
  }));
}

/**
 * Extract multiple ADMET results from structured_data.embedded_objects.
 */
function extractMultiAdmetFromStructured(response: AgentChatResponse): ADMETPropertyResult[] | undefined {
  const objects = response.structured_data?.embedded_objects;
  if (!objects?.length) return undefined;
  
  const admetObjects = objects.filter(obj => obj.type === 'admet_structured');
  if (!admetObjects.length) return undefined;
  
  return admetObjects.map(admetObj => ({
    smiles: admetObj.smiles || '',
    predictions: admetObj.predictions || {},
  }));
}

// =============================================================================
// UTILITY FUNCTIONS
// =============================================================================

/**
 * Check if a response contains structured data that needs visualization.
 * 
 * @param response - Raw response from backend
 * @returns true if response has structured data beyond text
 */
export function isStructuredResponse(response: AgentChatResponse): boolean {
  const responseType = response.response_type || 'text_only';
  return STRUCTURED_RESPONSE_TYPES.includes(responseType);
}

/**
 * Get a human-readable description of the response type.
 * Useful for UI labels and debugging.
 * 
 * @param responseType - The response type from backend
 * @returns Human-readable description
 */
export function getResponseTypeLabel(responseType: ChatResponseType): string {
  const labels: Record<ChatResponseType, string> = {
    'text_only':         'Text Response',
    'annotation':        'Molecule Analysis',
    'query_set':         'Query Set',
    'bioisosteres':      'Bioisostere Results',
    'admet':             'ADMET Predictions',
    'multi_bioisostere': 'Multiple Bioisostere Results',
    'multi_admet':       'Multiple ADMET Predictions',
  };
  return labels[responseType] || 'Unknown';
}

/**
 * Check if a parsed response has valid data for its action type.
 * Useful for defensive rendering.
 * 
 * @param parsed - The parsed response
 * @returns true if the expected data field is populated
 */
export function hasValidDataForAction(parsed: ParsedChatResponse): boolean {
  switch (parsed.action) {
    case 'RENDER_TEXT':
      return true; // Text is always valid
    case 'RENDER_MOLECULE':
      return !!parsed.moleculeData?.svg;
    case 'RENDER_QUERY_SET':
      return !!parsed.querySetData?.targets?.length;
    case 'RENDER_BIOISOSTERES':
      return !!parsed.bioisostereData?.results;
    case 'RENDER_ADMET':
      return !!parsed.admetData?.predictions;
    case 'RENDER_MULTI_BIO':
      return !!parsed.multiBioisostereData?.length;
    case 'RENDER_MULTI_ADMET':
      return !!parsed.multiAdmetData?.length;
    default:
      return false;
  }
}

// =============================================================================
// DEBUG UTILITIES
// =============================================================================

/**
 * Log a parsed response for debugging.
 * Color-coded console output similar to backend debug logging.
 */
export function debugLogParsedResponse(parsed: ParsedChatResponse): void {
  const colors = {
    action: 'color: #4CAF50; font-weight: bold',
    type: 'color: #2196F3',
    data: 'color: #FF9800',
    text: 'color: #9E9E9E',
  };
  
  console.group('%c📦 Parsed Chat Response', colors.action);
  console.log('%cAction:', colors.type, parsed.action);
  console.log('%cResponse Type:', colors.type, parsed.responseType);
  console.log('%cHas Structured Data:', colors.type, parsed.hasStructuredData);
  console.log('%cText Content:', colors.text, parsed.textContent.substring(0, 100) + '...');
  
  if (parsed.moleculeData) {
    console.log('%cMolecule Data:', colors.data, {
      smiles: parsed.moleculeData.smiles,
      matchCount: parsed.moleculeData.matches?.length || 0,
    });
  }
  if (parsed.bioisostereData) {
    console.log('%cBioisostere Data:', colors.data, {
      query: parsed.bioisostereData.query_smiles,
      resultCount: parsed.bioisostereData.num_results,
    });
  }
  if (parsed.querySetData) {
    console.log('%cQuery Set Data:', colors.data, {
      targetCount: parsed.querySetData.targets?.length || 0,
      queryType: parsed.querySetData.query_type,
    });
  }
  if (parsed.admetData) {
    console.log('%cADMET Data:', colors.data, {
      smiles: parsed.admetData.smiles,
      propertyCount: Object.keys(parsed.admetData.predictions || {}).length,
    });
  }
  
  console.groupEnd();
}

