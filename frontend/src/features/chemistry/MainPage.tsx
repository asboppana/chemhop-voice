import React, { useState, useEffect } from 'react';
import { useStructuredData, setMoleculeContext, clearMoleculeContext } from '@/hooks/useStructuredData';
import type { ParsedChatResponse } from '@/utils/chatResponseHandler';
import { FloatingNullState } from '@/components/animations/FloatingNullState';
import { analyzeMolecule, highlightMolecule, scanBioisosteres, generateSvg, llmReplaceSubstructure, predictADMET } from '@/services/moleculeService';
import type { MoleculeAnalysisResponse, BioisostereScanResponse, ADMETPredictionResponse } from '@/services/moleculeService';

// Helper function to convert index to letter (0 -> A, 1 -> B, etc.)
const indexToLetter = (index: number): string => {
  return String.fromCharCode(65 + index);
};

// Example molecules for quick access
const EXAMPLE_MOLECULES = [
  { name: 'Sildenafil', description: 'PDE5 inhibitor (Viagra)', smiles: 'CCCC1=NN(C2=C1N=C(NC2=O)C3=C(C=CC(=C3)S(=O)(=O)N4CCN(CC4)C)OCC)C' },
  { name: 'Ribociclib', description: 'CDK4/6 inhibitor', smiles: 'CN1CCN(CC1)c2nc(Nc3ccc(cc3)N4CCNCC4)nc(n2)C5CCCCC5' },
  { name: 'Lapatinib', description: 'Dual kinase inhibitor', smiles: 'CS(=O)(=O)CCNCc1ccc(cc1Cl)c2ccc3ncnc(Oc4ccc(F)c(Cl)c4)c3c2' },
  { name: 'Osimertinib', description: 'Multi-kinase inhibitor', smiles: 'CN1C=C(C2=CC=CC=C21)C3=NC(=NC=C3)NC4=C(C=C(C(=C4)NC(=O)C=C)N(C)CCN(C)C)OC' },
  { name: 'Osimertinib', description: 'EGFR inhibitor', smiles: 'COc1cc(N(C)CCN(C)C)c(NC(=O)C=C)cc1Nc2nccc(n2)c3cn(C4CC4)c4ccccc34' },
  { name: 'Erlotinib', description: 'EGFR inhibitor', smiles: 'COCCOc1cc2c(Nc3cccc(c3)C#C)ncnc2cc1OCCOC' },
];

export const ChemistryMainPage: React.FC = () => {
  const [analyzedSmiles, setAnalyzedSmiles] = useState(''); // Store the original SMILES used for analysis
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [result, setResult] = useState<MoleculeAnalysisResponse | null>(null);
  const [displayedSvg, setDisplayedSvg] = useState<string | null>(null);
  const [highlightLoading, setHighlightLoading] = useState(false);
  
  // Multi-selection and expression state
  const [selectedPatterns, setSelectedPatterns] = useState<Set<number>>(new Set());
  const [expressionGroups, setExpressionGroups] = useState<Array<{pattern: number, variants: number}>>([]); // Array of single pattern queries
  const [variantCount, setVariantCount] = useState<number>(10);
  
  // Bio-isostere scanning state
  const [bioisostereResults, setBioisostereResults] = useState<Map<number, BioisostereScanResponse>>(new Map());
  const [scanningQueries, setScanningQueries] = useState(false);
  const [bioisostereSvgs, setBioisostereSvgs] = useState<Map<string, string>>(new Map()); // SMILES -> SVG
  const [loadingSvgs, setLoadingSvgs] = useState(false);
  
  // LLM-generated molecule state
  const [generatedMolecule, setGeneratedMolecule] = useState<{
    smiles: string;
    svg: string;
    sourceFragment: string;
    replacementFragment: string;
    explanation: string | null;
  } | null>(null);
  const [generatingMolecule, setGeneratingMolecule] = useState(false);
  
  // ADMET predictions state
  const [admetPredictions, setAdmetPredictions] = useState<ADMETPredictionResponse | null>(null);
  const [loadingAdmet, setLoadingAdmet] = useState(false);

  // ===========================================================================
  // STRUCTURED DATA FROM CHAT
  // ===========================================================================
  // Subscribe to chat response events and handle them appropriately
  // Uses the onData callback pattern - latestData is available but not needed here
  useStructuredData({
    onData: handleStructuredDataFromChat,
  });

  /**
   * Handle structured data from chat responses.
   * Maps each response type to the appropriate state update.
   * 
   * Response Types → Actions:
   * - RENDER_MOLECULE     → setResult with molecule annotation
   * - RENDER_BIOISOSTERES → setBioisostereResults with scan results
   * - RENDER_QUERY_SET    → (optional) show pending query targets
   * - RENDER_ADMET        → (future) show ADMET predictions
   */
  function handleStructuredDataFromChat(data: ParsedChatResponse) {
    console.log('🧪 MainPage received structured data:', data.action);
    
    switch (data.action) {
      case 'RENDER_MOLECULE':
        if (data.moleculeData) {
          // Update molecule visualization
          setResult(data.moleculeData);
          setAnalyzedSmiles(data.moleculeData.smiles);
          setDisplayedSvg(data.moleculeData.svg);
          // Reset other state
          setSelectedPatterns(new Set());
          setExpressionGroups([]);
          setBioisostereResults(new Map());
          setError(null);
          console.log('✅ Molecule rendered:', data.moleculeData.smiles);
        }
        break;
        
      case 'RENDER_BIOISOSTERES':
        if (data.bioisostereData) {
          // Convert to the Map format expected by the UI
          // For chat responses, we use index 0 as the key
          const newResults = new Map<number, BioisostereScanResponse>();
          newResults.set(0, {
            query_smiles: data.bioisostereData.query_smiles,
            num_matches: data.bioisostereData.num_results,
            matches: data.bioisostereData.results.map(r => ({
              source: r.source,
              similarity: r.similarity,
              centroid_smiles: r.centroid_smiles,
              bio_isostere_score: r.bio_isostere_score,
              pharmacophore_similarity: r.pharmacophore_similarity || 0,
              topology_similarity: r.topology_similarity || 0,
              cluster_id: undefined,
              num_members: 1,
              example_smiles: [],
              descriptors: r.descriptors,
              delta_properties: r.delta_properties || {},
            })),
          });
          setBioisostereResults(newResults);
          console.log('✅ Bioisosteres rendered:', data.bioisostereData.num_results, 'results');
        }
        break;
        
      case 'RENDER_MULTI_BIO':
        if (data.multiBioisostereData) {
          // Convert array of results to Map
          const newResults = new Map<number, BioisostereScanResponse>();
          data.multiBioisostereData.forEach((scanResult, index) => {
            newResults.set(index, {
              query_smiles: scanResult.query_smiles,
              num_matches: scanResult.num_results,
              matches: scanResult.results.map(r => ({
                source: r.source,
                similarity: r.similarity,
                centroid_smiles: r.centroid_smiles,
                bio_isostere_score: r.bio_isostere_score,
                pharmacophore_similarity: r.pharmacophore_similarity || 0,
                topology_similarity: r.topology_similarity || 0,
                cluster_id: undefined,
                num_members: 1,
                example_smiles: [],
                descriptors: r.descriptors,
                delta_properties: r.delta_properties || {},
              })),
            });
          });
          setBioisostereResults(newResults);
          console.log('✅ Multi-bioisosteres rendered:', newResults.size, 'query results');
        }
        break;
        
      case 'RENDER_QUERY_SET':
        if (data.querySetData) {
          // Query set shows targets that will be queried
          // Could display as pending state or just log
          console.log('📋 Query set received:', data.querySetData.targets?.length, 'targets');
          // Optionally: set some UI state to show pending queries
        }
        break;
        
      case 'RENDER_ADMET':
        if (data.admetData) {
          // TODO: Implement ADMET visualization
          console.log('📊 ADMET data received:', Object.keys(data.admetData.predictions).length, 'properties');
        }
        break;
        
      case 'RENDER_MULTI_ADMET':
        if (data.multiAdmetData) {
          // TODO: Implement multi-ADMET visualization
          console.log('📊 Multi-ADMET data received:', data.multiAdmetData.length, 'molecules');
        }
        break;
        
      case 'RENDER_TEXT':
      default:
        // Text-only responses don't affect MainPage visualization
        break;
    }
  }

  // Generate SVGs for bio-isostere results
  useEffect(() => {
    const generateBioisostereSvgs = async () => {
      if (bioisostereResults.size === 0) return;
      
      setLoadingSvgs(true);
      const newSvgs = new Map<string, string>();
      
      try {
        // Collect all unique SMILES
        const allSmiles = new Set<string>();
        bioisostereResults.forEach(scanResult => {
          scanResult.matches.forEach(match => {
            allSmiles.add(match.centroid_smiles);
          });
        });
        
        // Generate SVG for each unique SMILES
        for (const smiles of allSmiles) {
          try {
            const svgResponse = await generateSvg(smiles, 120, 120);
            newSvgs.set(smiles, svgResponse.svg);
          } catch (err) {
            console.error(`Error generating SVG for ${smiles}:`, err);
          }
        }
        
        setBioisostereSvgs(newSvgs);
      } catch (err) {
        console.error('Error generating bio-isostere SVGs:', err);
      } finally {
        setLoadingSvgs(false);
      }
    };
    
    generateBioisostereSvgs();
  }, [bioisostereResults]);

  // ===========================================================================
  // PUBLISH MOLECULE CONTEXT FOR CHAT
  // ===========================================================================
  // When molecule analysis changes, publish the context so ChatContext can
  // include it in conversation history for the backend to reference
  useEffect(() => {
    if (!result || !analyzedSmiles) {
      clearMoleculeContext();
      return;
    }

    // Filter for cyclic/biological/ring patterns (same filter used in UI)
    const filteredMatches = result.matches.filter(match => 
      match.trivial_name.group.toLowerCase().includes('cyclic') ||
      match.trivial_name.group.toLowerCase().includes('biological') ||
      match.trivial_name.group.toLowerCase().includes('ring')
    );

    // Build pattern info array
    const patterns = filteredMatches.map((match, index) => ({
      letter: indexToLetter(index),
      name: match.trivial_name.name,
      smiles: match.trivial_name.smarts,
      atomIndices: match.atom_indices,
    }));

    // Publish context
    setMoleculeContext({
      smiles: analyzedSmiles,
      name: result.name,
      patterns,
      timestamp: Date.now(),
    });

    // Cleanup when component unmounts
    return () => {
      clearMoleculeContext();
    };
  }, [result, analyzedSmiles]);

  const handleAnalyzeMolecule = async (smiles: string, moleculeName?: string) => {
    if (!smiles.trim()) {
      setError('Please enter a SMILES string');
      return;
    }

    setLoading(true);
    setError(null);
    
    try {
      const response = await analyzeMolecule(smiles);
      // If we have a provided name (e.g., from example molecules), use it
      // Otherwise use the backend's response name
      const resultWithName = moleculeName 
        ? { ...response, name: moleculeName }
        : response;
      setResult(resultWithName);
      setAnalyzedSmiles(smiles); // Store the original SMILES for highlighting
      setDisplayedSvg(response.svg);
      // Reset selection state
      setSelectedPatterns(new Set());
      setExpressionGroups([]);
      setBioisostereResults(new Map());
      console.log('Molecule analysis result:', resultWithName);
    } catch (err: any) {
      setError(err.response?.data?.detail || 'Failed to analyze molecule');
      console.error('Error analyzing molecule:', err);
    } finally {
      setLoading(false);
    }
  };

  // Toggle pattern selection (single select only)
  const handlePatternClick = async (_atomIndices: number[], patternIndex: number) => {
    if (!result) return;
    
    const newSelected = new Set<number>();
    if (!selectedPatterns.has(patternIndex)) {
      newSelected.add(patternIndex);
    }
    setSelectedPatterns(newSelected);
    
    // Highlight selected pattern
    if (newSelected.size > 0) {
      // Get filtered matches to access atom indices
      const filteredMatches = result.matches.filter(match => 
        match.trivial_name.group.toLowerCase().includes('cyclic') ||
        match.trivial_name.group.toLowerCase().includes('biological') ||
        match.trivial_name.group.toLowerCase().includes('ring')
      );
      
      const patternIdx = Array.from(newSelected)[0];
      if (filteredMatches[patternIdx]) {
        setHighlightLoading(true);
        try {
          const response = await highlightMolecule(analyzedSmiles, filteredMatches[patternIdx].atom_indices);
          setDisplayedSvg(response.svg);
        } catch (err: any) {
          console.error('Error highlighting molecule:', err);
        } finally {
          setHighlightLoading(false);
        }
      }
    } else {
      setDisplayedSvg(result.svg);
    }
  };

  // Handle + button - add single pattern as new query
  const handleAddQuery = () => {
    if (selectedPatterns.size === 0) return;
    
    const patternIdx = Array.from(selectedPatterns)[0]; // Only one pattern since we're single-select
    
    // Check if this pattern already exists in queries
    const isDuplicate = expressionGroups.some(group => group.pattern === patternIdx);
    
    if (isDuplicate) {
      console.log('Pattern already in queries, skipping');
      setSelectedPatterns(new Set());
      return;
    }
    
    // Add as a new query
    setExpressionGroups([...expressionGroups, { pattern: patternIdx, variants: variantCount }]);
    
    // Clear selection
    setSelectedPatterns(new Set());
  };

  // Build expression lines for display
  const buildExpressionLines = (): Array<{expression: string, variants: number}> => {
    if (expressionGroups.length === 0) return [];
    
    return expressionGroups.map(group => ({
      expression: indexToLetter(group.pattern),
      variants: group.variants
    }));
  };

  // Clear expression
  const handleClearExpression = () => {
    setExpressionGroups([]);
    setSelectedPatterns(new Set());
    setBioisostereResults(new Map());
    setBioisostereSvgs(new Map());
    setGeneratedMolecule(null);
    if (result) {
      setDisplayedSvg(result.svg);
    }
  };
  
  // Handle bio-isostere match click - use LLM to replace substructure
  const handleBioisostereClick = async (
    patternIdx: number,
    matchSmiles: string,
    _queryGroupLetter: string
  ) => {
    if (!result) return;
    
    // Get filtered matches to access the source pattern
    const filteredMatches = result.matches.filter(match => 
      match.trivial_name.group.toLowerCase().includes('cyclic') ||
      match.trivial_name.group.toLowerCase().includes('biological') ||
      match.trivial_name.group.toLowerCase().includes('ring')
    );
    
    const sourcePattern = filteredMatches[patternIdx];
    if (!sourcePattern) return;
    
    setGeneratingMolecule(true);
    setGeneratedMolecule(null);
    setAdmetPredictions(null);
    
    try {
      console.log('Calling LLM replace substructure...');
      console.log('Original:', analyzedSmiles);
      console.log('Source fragment:', sourcePattern.trivial_name.smarts);
      console.log('Replacement fragment:', matchSmiles);
      
      // Call the LLM replace substructure endpoint
      const replaceResponse = await llmReplaceSubstructure(
        analyzedSmiles,
        sourcePattern.trivial_name.smarts,
        matchSmiles
      );
      
      console.log('LLM Replace Response:', replaceResponse);
      
      if (replaceResponse.success && replaceResponse.result_smiles) {
        // Generate SVG for the result molecule
        const svgResponse = await generateSvg(replaceResponse.result_smiles, 400, 400);
        
        setGeneratedMolecule({
          smiles: replaceResponse.result_smiles,
          svg: svgResponse.svg,
          sourceFragment: sourcePattern.trivial_name.smarts,
          replacementFragment: matchSmiles,
          explanation: replaceResponse.explanation
        });
        
        // Call ADMET prediction for the generated molecule
        setLoadingAdmet(true);
        try {
          console.log('Calling ADMET prediction...');
          const admetResponse = await predictADMET(replaceResponse.result_smiles);
          console.log('ADMET Response:', admetResponse);
          setAdmetPredictions(admetResponse);
        } catch (admetErr: any) {
          console.error('Error getting ADMET predictions:', admetErr);
          // Don't set error - ADMET is optional, molecule generation succeeded
        } finally {
          setLoadingAdmet(false);
        }
      } else {
        console.warn('LLM replacement was not successful:', replaceResponse.explanation);
        setError(replaceResponse.explanation || 'Failed to generate replacement molecule');
      }
      
    } catch (err: any) {
      console.error('Error in LLM replacement:', err);
      setError(err.message || 'Failed to generate replacement molecule');
    } finally {
      setGeneratingMolecule(false);
    }
  };

  // Run bio-isostere queries for all patterns
  const handleRunQueries = async () => {
    if (!result || expressionGroups.length === 0) return;
    
    setScanningQueries(true);
    const newResults = new Map<number, BioisostereScanResponse>();
    
    try {
      // Get filtered matches (same filter as display)
      const filteredMatches = result.matches.filter(match => 
        match.trivial_name.group.toLowerCase().includes('cyclic') ||
        match.trivial_name.group.toLowerCase().includes('biological') ||
        match.trivial_name.group.toLowerCase().includes('ring')
      );
      
      // Run bio-isostere scan for each query
      for (const group of expressionGroups) {
        const patternIdx = group.pattern;
        const variants = group.variants;
        
        try {
          const match = filteredMatches[patternIdx];
          if (!match) continue;
          
          // Use the SMARTS/SMILES from the pattern directly
          // For cyclic patterns, the smarts field contains the SMILES string
          const ringSmiles = match.trivial_name.smarts;
          
          // Scan for bio-isosteres
          const scanResponse = await scanBioisosteres(
            ringSmiles,
            variants,
            0.3 // default min_similarity
          );
          
          newResults.set(patternIdx, scanResponse);
        } catch (err) {
          console.error(`Error scanning pattern ${indexToLetter(patternIdx)}:`, err);
        }
      }
      
      setBioisostereResults(newResults);
    } catch (err: any) {
      console.error('Error running queries:', err);
      setError(err.response?.data?.detail || 'Failed to run bio-isostere queries');
    } finally {
      setScanningQueries(false);
    }
  };

  return (
    <div className="w-full bg-white text-black">
      {/* Main content area */}
      <div className="min-h-full">
        {loading && (
          <div className="flex flex-col items-center justify-center min-h-[60vh] py-20">
            <FloatingNullState compact={true} />
            <p className="text-sm text-gray-600 mt-8">Analyzing molecule...</p>
          </div>
        )}

        {error && (
          <div className="flex items-center justify-center min-h-[60vh] py-20">
            <div className="text-center max-w-md">
              <p className="text-sm text-red-600 mb-2">Error</p>
              <p className="text-xs text-gray-600">{error}</p>
            </div>
          </div>
        )}

        {!loading && !error && !result && (
          <div className="flex flex-col items-center justify-start w-full pt-16 pb-8">
            {/* Welcome header */}
            <div className="text-center mb-10">
              <div className="w-16 h-16 mx-auto mb-5 rounded-xl bg-gray-100 flex items-center justify-center border border-gray-500">
                <img src="/icons/Beaker.svg" alt="Beaker" className="w-8 h-8" />
              </div>
              <h1 className="text-2xl font-medium text-black mb-3">
                ChemHop Voice
              </h1>
              <p className="text-base text-gray-1000">
                Analyze molecular structures and discover bioisostere replacements
              </p>
            </div>

            {/* Example molecules grid */}
            <div className="w-full max-w-lg px-4">
              <p className="text-xs uppercase tracking-[0.02em] text-gray-600 mb-4 text-center font-medium">
                Try an example
              </p>
              <div className="grid grid-cols-2 gap-3">
                {EXAMPLE_MOLECULES.map((example) => (
                  <button
                    key={example.smiles}
                    onClick={() => handleAnalyzeMolecule(example.smiles, example.name)}
                    className="group text-left px-4 py-3.5 rounded-lg border border-gray-200 hover:border-gray-300 hover:scale-[0.98] hover:shadow-md transition-all duration-300"
                  >
                    <p className="text-sm font-semibold text-gray-900 group-hover:text-black">
                      {example.name}
                    </p>
                    <p className="text-sm text-gray-600 mt-1">
                      {example.description}
                    </p>
                  </button>
                ))}
              </div>
            </div>            
          </div>
        )}

        {!loading && !error && result && (
          <div className="space-y-6 p-12 pb-32">
            {/* Source Molecule and Patterns Side by Side */}
            <div className="flex flex-col lg:flex-row gap-6 items-start">
              {/* SVG Visualization - Larger box on left */}
              <div className="flex-shrink-0 w-full lg:w-auto">
                <div className="flex items-center h-8 mb-3">
                  <h3 className="text-sm font-medium uppercase tracking-[0.02em] text-black">
                    {result.name && result.name !== 'No Name' && result.name !== 'Unknown' ? (
                      <>{result.name} <span className="text-[9px] pl-1 text-gray-1000">[Source Molecule]</span></>
                    ) : (
                      'Source Molecule'
                    )}
                  </h3>
                </div>
                <div className="border border-gray-300 rounded-lg p-4 bg-white w-full max-w-md mx-auto lg:mx-0 aspect-square lg:w-96 relative">
                  {highlightLoading && (
                    <div className="absolute inset-0 bg-white/70 flex items-center justify-center z-10">
                      <div className="animate-spin rounded-full h-8 w-8 border-b-2 border-black"></div>
                    </div>
                  )}
                  <div 
                    dangerouslySetInnerHTML={{ __html: displayedSvg || result.svg }}
                    className="pattern-svg-container w-full h-full flex items-center justify-center"
                  />
                </div>
                <p className="text-[10px] text-gray-600 font-mono mt-2 break-all max-w-md mx-auto lg:mx-0">{result.smiles}</p>
                
                {/* Query display */}
                {expressionGroups.length > 0 && (
                  <div className="mt-4 p-3 bg-gray-100 rounded-lg max-w-md mx-auto lg:mx-0">
                    <p className="text-[9px] uppercase tracking-[0.02em] text-gray-1500 font-medium mb-1.5">Queries</p>
                    <div className="space-y-0.5 mb-3">
                      {buildExpressionLines().map((line, idx) => (
                        <div key={idx} className="flex items-baseline gap-1.5">
                          <span className="text-sm font-medium text-gray-1000 w-6 text-left flex-shrink-0">{idx + 1}.</span>
                          <span className="text-sm font-mono font-medium text-black flex-1">{line.expression}</span>
                          <span className="text-[10px] text-gray-600 font-mono flex-shrink-0">[{line.variants}]</span>
                        </div>
                      ))}
                    </div>
                    <button
                      onClick={handleRunQueries}
                      disabled={scanningQueries}
                      className={`w-full py-2 text-[11px] font-semibold uppercase tracking-wider rounded transition-all flex items-center justify-center gap-2 ${
                        scanningQueries
                          ? 'bg-gray-200 cursor-not-allowed'
                          : 'bg-black text-white hover:bg-gray-800'
                      }`}
                    >
                      {scanningQueries && (
                        <div className="animate-spin rounded-full h-3 w-3 border-2 border-gray-1500 border-t-transparent"></div>
                      )}
                      <span className={scanningQueries ? 'text-gray-1500' : ''}>
                        {scanningQueries ? 'Scanning...' : 'Run Queries'}
                      </span>
                    </button>
                  </div>
                )}
              </div>

              {/* Matched Patterns as mini boxes - filter for cyclic/biological */}
              {result.matches && result.matches.length > 0 && (() => {
                const filteredMatches = result.matches.filter(match => 
                  match.trivial_name.group.toLowerCase().includes('cyclic') ||
                  match.trivial_name.group.toLowerCase().includes('biological') ||
                  match.trivial_name.group.toLowerCase().includes('ring')
                );
                
                return filteredMatches.length > 0 ? (
                  <div className="flex-1 min-w-0">
                    {/* Header with + button */}
                    <div className="flex items-center justify-between h-8 mb-3">
                      <h3 className="text-sm font-medium uppercase tracking-[0.02em] text-black">
                        Patterns <span className="text-sm pl-1 text-gray-1000">[{filteredMatches.length}]</span>
                      </h3>
                      <div className="flex items-center gap-2.5">
                        {/* Variant count input */}
                        <div className="flex items-center gap-2 border border-gray-200 rounded-lg px-3 py-1.5 bg-white">
                          <label htmlFor="variant-count" className="text-xs uppercase tracking-[0.02em] text-gray-1500 font-medium">
                            Variants
                          </label>
                          <input
                            id="variant-count"
                            type="number"
                            min="1"
                            max="100"
                            value={variantCount}
                            onChange={(e) => setVariantCount(Math.max(1, Math.min(100, parseInt(e.target.value) || 1)))}
                            className="w-10 bg-transparent text-sm font-medium text-black text-center focus:outline-none"
                          />
                        </div>
                        
                        <button
                          onClick={handleAddQuery}
                          disabled={selectedPatterns.size === 0}
                          className={`w-7 h-7 rounded-full font-semibold text-base tracking-[0.02em] transition-all flex items-center justify-center ${
                            selectedPatterns.size > 0
                              ? 'bg-black text-white hover:bg-gray-800'
                              : 'bg-gray-200 text-gray-500 cursor-not-allowed'
                          }`}
                          title="Add pattern to queries"
                        >
                          +
                        </button>
                        {expressionGroups.length > 0 && (
                          <button
                            onClick={handleClearExpression}
                            className="w-7 h-7 rounded-full bg-black hover:bg-gray-800 text-sm text-white flex items-center justify-center transition-colors"
                            title="Clear queries"
                          >
                            ✕
                          </button>
                        )}
                      </div>
                    </div>
                    
                    {/* Pattern grid */}
                    <div className="grid gap-2" style={{ gridTemplateColumns: 'repeat(auto-fill, minmax(140px, 1fr))' }}>
                      {filteredMatches.map((match, index) => (
                        <button 
                          key={index}
                          onClick={() => handlePatternClick(match.atom_indices, index)}
                          className={`border rounded p-2 transition-colors flex flex-col aspect-square cursor-pointer text-left relative w-full ${
                            selectedPatterns.has(index) 
                              ? 'border-blue-500 bg-blue-50 ring-2 ring-blue-200' 
                              : 'border-gray-300 bg-gray-50 hover:border-gray-400'
                          }`}
                        >
                          {/* Letter label badge */}
                          <div className="absolute top-1 left-1 w-5 h-5 rounded-full bg-black text-white text-[10px] font-medium flex items-center justify-center">
                            {indexToLetter(index)}
                          </div>
                          <p className="text-xs font-medium mb-1 pl-6 leading-tight min-h-[28px]">{match.trivial_name.name}</p>
                          <div className="flex-1 min-h-0 flex items-center justify-center p-1">
                            {match.svg ? (
                              <div 
                                dangerouslySetInnerHTML={{ __html: match.svg }}
                                className="pattern-svg-container w-full h-full flex items-center justify-center"
                              />
                            ) : (
                              <p className="text-[8px] text-gray-400 text-center break-all leading-tight p-2">
                                {match.trivial_name.smarts}
                              </p>
                            )}
                          </div>
                        </button>
                      ))}
                    </div>
                  </div>
                ) : null;
              })()}
            </div>

            {/* Generated Molecule Result */}
            {generatedMolecule && (
              <div className="mt-12 border-t border-gray-300 pt-6">
                <h2 className="text-sm font-medium mb-6 uppercase tracking-wider text-black">
                  Generated Molecule
                </h2>
                
                <div className="flex flex-col lg:flex-row gap-6">
                  {/* Molecule View - Left Side */}
                  <div className="flex-shrink-0 w-full lg:w-auto">
                    <div className="border border-gray-300 rounded-lg p-4 bg-white w-full max-w-md mx-auto lg:mx-0 aspect-square lg:w-96">
                      <div 
                        dangerouslySetInnerHTML={{ __html: generatedMolecule.svg }}
                        className="pattern-svg-container w-full h-full flex items-center justify-center"
                      />
                    </div>
                    <p className="text-[10px] text-gray-600 font-mono mt-2 break-all max-w-md mx-auto lg:mx-0">
                      {generatedMolecule.smiles}
                    </p>
                    
                    {/* Replacement Info */}
                    <div className="mt-4 p-3 bg-gray-50 rounded-lg max-w-md mx-auto lg:mx-0">
                      <p className="text-[9px] uppercase tracking-wider text-gray-600 mb-2">Replacement Details</p>
                      <div className="space-y-1.5">
                        <div className="flex items-start gap-2">
                          <span className="text-[10px] text-gray-500 w-16 flex-shrink-0">Source:</span>
                          <span className="text-[10px] font-mono text-gray-700 break-all">{generatedMolecule.sourceFragment}</span>
                        </div>
                        <div className="flex items-start gap-2">
                          <span className="text-[10px] text-gray-500 w-16 flex-shrink-0">Replaced:</span>
                          <span className="text-[10px] font-mono text-gray-700 break-all">{generatedMolecule.replacementFragment}</span>
                        </div>
                      </div>
                      {generatedMolecule.explanation && (
                        <p className="text-[10px] text-gray-600 mt-3 pt-2 border-t border-gray-200">
                          {generatedMolecule.explanation}
                        </p>
                      )}
                    </div>
                  </div>
                  
                  {/* ADMET Properties - Right Side */}
                  <div className="flex-1 min-w-0">
                    <h3 className="text-sm font-medium uppercase tracking-wider text-black mb-3">
                      ADMET Properties
                    </h3>
                    <div className="border border-gray-300 rounded-lg bg-white overflow-hidden">
                      {loadingAdmet ? (
                        <div className="flex items-center justify-center py-12">
                          <div className="animate-spin rounded-full h-8 w-8 border-b-2 border-gray-600"></div>
                          <span className="ml-3 text-sm text-gray-600">Loading ADMET predictions...</span>
                        </div>
                      ) : admetPredictions?.error ? (
                        <div className="px-4 py-8 text-center">
                          <p className="text-sm text-gray-500">{admetPredictions.error}</p>
                        </div>
                      ) : admetPredictions?.predictions && Object.keys(admetPredictions.predictions).length > 0 ? (
                        <div className="max-h-[500px] overflow-y-auto">
                          <table className="w-full text-sm">
                            <thead className="bg-gray-50 border-b border-gray-200 sticky top-0">
                              <tr>
                                <th className="px-4 py-2 text-left text-xs font-medium text-gray-600 uppercase tracking-wider">Property</th>
                                <th className="px-4 py-2 text-right text-xs font-medium text-gray-600 uppercase tracking-wider">Value</th>
                              </tr>
                            </thead>
                            <tbody className="divide-y divide-gray-100">
                              {Object.entries(admetPredictions.predictions)
                                .sort(([a], [b]) => a.localeCompare(b))
                                .map(([property, value]) => (
                                  <tr key={property} className="hover:bg-gray-50">
                                    <td className="px-4 py-2 text-gray-700 text-xs">
                                      {property.replace(/_/g, ' ')}
                                    </td>
                                    <td className="px-4 py-2 text-right font-mono text-gray-600 text-xs">
                                      {typeof value === 'number' ? value.toFixed(3) : String(value)}
                                    </td>
                                  </tr>
                                ))
                              }
                            </tbody>
                          </table>
                        </div>
                      ) : (
                        <div className="px-4 py-8 text-center">
                          <p className="text-sm text-gray-500">No ADMET predictions available</p>
                          <p className="text-xs text-gray-400 mt-1">ADMET service may not be running</p>
                        </div>
                      )}
                    </div>
                  </div>
                </div>
              </div>
            )}

            {/* Bio-isostere Results */}
            {bioisostereResults.size > 0 && result && (() => {
              const filteredMatches = result.matches.filter(match => 
                match.trivial_name.group.toLowerCase().includes('cyclic') ||
                match.trivial_name.group.toLowerCase().includes('biological') ||
                match.trivial_name.group.toLowerCase().includes('ring')
              );

              return (
                <div className="mt-12">
                  <h2 className="text-lg font-medium mb-6 uppercase tracking-[0.02em] text-black">
                    Bioisostere Results
                  </h2>
                  
                  <div className="space-y-8">
                    {Array.from(bioisostereResults.entries()).map(([patternIdx, scanResult]) => {
                      const originalPattern = filteredMatches[patternIdx];
                      if (!originalPattern) return null;
                      
                      return (
                        <div key={patternIdx} className="border-t border-gray-300 pt-6">
                          <div className="mb-4">
                            {/* Source Pattern Display */}
                            <div className="inline-block border border-gray-300 rounded p-2 bg-gray-50 w-[140px]">
                              <div className="flex items-start gap-1.5 mb-2 min-h-[32px]">
                                <div className="w-5 h-5 rounded-full bg-black text-white text-[10px] font-medium flex items-center justify-center flex-shrink-0">
                                  {indexToLetter(patternIdx)}
                                </div>
                                <span className="text-[10px] font-medium text-black leading-tight">{originalPattern.trivial_name.name}</span>
                              </div>
                              <div className="w-full h-12 flex items-center justify-center mb-2 overflow-hidden">
                                {originalPattern.svg ? (
                                  <div 
                                    dangerouslySetInnerHTML={{ __html: originalPattern.svg }}
                                    className="w-full h-full [&>svg]:w-full [&>svg]:h-full [&>svg]:object-contain"
                                    style={{ display: 'flex', alignItems: 'center', justifyContent: 'center' }}
                                  />
                                ) : (
                                  <p className="text-[8px] text-gray-600 font-mono text-center break-all leading-tight">
                                    {originalPattern.trivial_name.smarts}
                                  </p>
                                )}
                              </div>
                              <p className="text-[9px] text-gray-1500 font-mono truncate">
                                {scanResult.query_smiles}
                              </p>
                              <p className="text-[9px] text-gray-600 mt-1">
                                {scanResult.num_matches} matches
                              </p>
                            </div>
                          </div>
                          
                          {/* Matches grid */}
                          <div className="grid gap-2" style={{ gridTemplateColumns: 'repeat(auto-fill, minmax(120px, 1fr))' }}>
                            {scanResult.matches.map((match, matchIdx) => {
                              const svg = bioisostereSvgs.get(match.centroid_smiles);
                              const queryGroupLetter = indexToLetter(patternIdx);
                              
                              return (
                                <button
                                  key={matchIdx}
                                  onClick={() => handleBioisostereClick(patternIdx, match.centroid_smiles, queryGroupLetter)}
                                  disabled={generatingMolecule}
                                  className={`border border-gray-300 rounded p-2 bg-white transition-colors flex flex-col w-full text-left ${
                                    generatingMolecule 
                                      ? 'cursor-not-allowed opacity-50' 
                                      : 'hover:border-gray-400 hover:shadow-md cursor-pointer'
                                  }`}
                                  title={`Click to generate molecule with this bio-isostere replacement`}
                                >
                                  {/* Header */}
                                  <div className="flex items-center justify-between mb-1">
                                    <span className="text-[8px] font-semibold text-gray-1500 uppercase tracking-wider leading-none">
                                      Weighted<br/>Score
                                    </span>
                                    <span className="text-sm font-mono text-black font-medium">
                                      {match.bio_isostere_score.toFixed(2)}
                                    </span>
                                  </div>
                                  
                                  {/* SVG visualization */}
                                  <div className="w-full h-14 flex items-center justify-center mb-1 overflow-hidden">
                                    {loadingSvgs ? (
                                      <div className="animate-spin rounded-full h-6 w-6 border-b border-gray-400"></div>
                                    ) : svg ? (
                                      <div 
                                        dangerouslySetInnerHTML={{ __html: svg }}
                                        className="w-full h-full [&>svg]:w-full [&>svg]:h-full [&>svg]:object-contain"
                                        style={{ display: 'flex', alignItems: 'center', justifyContent: 'center' }}
                                      />
                                    ) : (
                                      <p className="text-[8px] text-gray-1500 text-center break-all leading-tight px-1">
                                        {match.centroid_smiles}
                                      </p>
                                    )}
                                  </div>
                                  
                                  {/* Footer - scores */}
                                  <div className="pt-1 border-t border-gray-200 space-y-0">
                                    <div className="flex justify-between text-[8px] leading-relaxed">
                                      <span className="text-gray-1500">Morgan:</span>
                                      <span className="font-mono text-gray-600">{match.similarity.toFixed(2)}</span>
                                    </div>
                                    <div className="flex justify-between text-[8px] leading-relaxed">
                                      <span className="text-gray-1500">Pharm:</span>
                                      <span className="font-mono text-gray-600">{match.pharmacophore_similarity.toFixed(2)}</span>
                                    </div>
                                    <div className="flex justify-between text-[8px] leading-relaxed">
                                      <span className="text-gray-1500">Topo:</span>
                                      <span className="font-mono text-gray-600">{match.topology_similarity.toFixed(2)}</span>
                                    </div>
                                    <div className="flex justify-between text-[8px] leading-relaxed pt-0.5 border-t border-gray-100">
                                      <span className="text-gray-1500">Source:</span>
                                      <span className="font-mono text-gray-600 capitalize">{match.source}</span>
                                    </div>
                                  </div>
                                </button>
                              );
                            })}
                          </div>
                        </div>
                      );
                    })}
                  </div>
                </div>
              );
            })()}
          </div>
        )}
      </div>
    </div>
  );
};

export default ChemistryMainPage;
