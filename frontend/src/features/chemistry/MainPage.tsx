import React, { useState, useEffect } from 'react';
import { analyzeMolecule, highlightMolecule, scanBioisosteres, generateSvg } from '@/services/moleculeService';
import type { MoleculeAnalysisResponse, BioisostereScanResponse } from '@/services/moleculeService';

// Helper function to convert index to letter (0 -> A, 1 -> B, etc.)
const indexToLetter = (index: number): string => {
  return String.fromCharCode(65 + index);
};

export const ChemistryMainPage: React.FC = () => {
  const [smilesInput, setSmilesInput] = useState('');
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

  const handleSubmit = async (e: React.FormEvent) => {
    e.preventDefault();
    
    if (!smilesInput.trim()) {
      setError('Please enter a SMILES string');
      return;
    }

    setLoading(true);
    setError(null);
    
    try {
      const response = await analyzeMolecule(smilesInput);
      setResult(response);
      setAnalyzedSmiles(smilesInput); // Store the original SMILES for highlighting
      setDisplayedSvg(response.svg);
      // Reset selection state
      setSelectedPatterns(new Set());
      setExpressionGroups([]);
      console.log('Molecule analysis result:', response);
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
    if (result) {
      setDisplayedSvg(result.svg);
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
    <div className="h-screen w-full bg-white text-black flex">
      {/* Left side - SMILES Input */}
      <div className="w-80 border-r border-gray-300 flex flex-col overflow-hidden">
        <div className="flex-shrink-0 p-6 pb-0">
          <h1 className="text-xl font-light mb-6 tracking-wide">SMILES INPUT</h1>
          
          <form onSubmit={handleSubmit} className="flex flex-col gap-4">
            <div className="flex flex-col gap-2">
              <label htmlFor="smiles" className="text-xs uppercase tracking-wider text-gray-600">
                Enter SMILES String
              </label>
              <input
                id="smiles"
                type="text"
                value={smilesInput}
                onChange={(e) => setSmilesInput(e.target.value)}
                placeholder="CCO"
                className="bg-transparent border border-gray-300 rounded px-3 py-2 text-sm focus:outline-none focus:border-black transition-colors font-mono"
              />
            </div>
            
            <button
              type="submit"
              className="bg-black text-white px-4 py-2 text-sm font-medium hover:bg-gray-800 transition-colors"
            >
              SUBMIT
            </button>
          </form>
        </div>

        {/* Example SMILES */}
        <div className="flex-1 min-h-0 mt-8 px-6 pb-6">
          <div className="border-t border-gray-300 pt-6">
            <p className="text-xs uppercase tracking-wider text-gray-600 mb-3">Examples</p>
            <div className="flex flex-col gap-2 overflow-y-auto" style={{ maxHeight: 'calc(100vh - 400px)' }}>
              {[
                { name: 'Sildenafil (Viagra)', smiles: 'CCCC1=NN(C2=C1N=C(NC2=O)C3=C(C=CC(=C3)S(=O)(=O)N4CCN(CC4)C)OCC)C' },
                { name: 'Ribociclib', smiles: 'CN1CCN(CC1)c2nc(Nc3ccc(cc3)N4CCNCC4)nc(n2)C5CCCCC5' },
                { name: 'Lapatinib', smiles: 'CS(=O)(=O)CCNCc1ccc(cc1Cl)c2ccc3ncnc(Oc4ccc(F)c(Cl)c4)c3c2' },
                { name: 'Sunitinib', smiles: 'CCN(CC)CCNC(=O)c1c(C)[nH]c(c1C)c2c[nH]c(=O)c(c2)c3cc(F)ccc3' },
                { name: 'Osimertinib', smiles: 'COc1cc(N(C)CCN(C)C)c(NC(=O)C=C)cc1Nc2nccc(n2)c3cn(C4CC4)c4ccccc34' },
                { name: 'Erlotinib', smiles: 'COCCOc1cc2c(Nc3cccc(c3)C#C)ncnc2cc1OCCOC' },
                { name: 'Palbociclib', smiles: 'CC(C)n1cc(cn1)c2cnc(nc2N3CCNCC3)Nc4ccc(cn4)C(=O)NC5CCCCC5' },
                { name: 'Venetoclax', smiles: 'CC1(C)CCC(C)(C)c2c1ccc(c2)c3cc(c(N4CCN(CC4)CCO)c(c3)c5ccc(cc5)c6cccnc6)C(=O)NS(=O)(=O)c7ccc(cc7)NCC8(CC8)C(F)(F)F' },
              ].map((example) => (
                <button
                  key={example.smiles}
                  onClick={() => setSmilesInput(example.smiles)}
                  className="text-left p-2 hover:bg-gray-100 transition-colors rounded text-xs flex-shrink-0"
                >
                  <div className="font-medium">{example.name}</div>
                  <div className="text-gray-600 font-mono text-[10px] mt-1">{example.smiles}</div>
                </button>
              ))}
            </div>
          </div>
        </div>
      </div>

      {/* Right side - Main content area */}
      <div className="flex-1 p-12 overflow-auto">
        {loading && (
          <div className="flex items-center justify-center h-full">
            <div className="text-center">
              <div className="animate-spin rounded-full h-12 w-12 border-b-2 border-black mx-auto mb-4"></div>
              <p className="text-sm text-gray-600">Analyzing molecule...</p>
            </div>
          </div>
        )}

        {error && (
          <div className="flex items-center justify-center h-full">
            <div className="text-center max-w-md">
              <p className="text-sm text-red-600 mb-2">Error</p>
              <p className="text-xs text-gray-600">{error}</p>
            </div>
          </div>
        )}

        {!loading && !error && !result && (
          <div className="flex items-center justify-center h-full">
            <div className="text-center">
              <p className="text-sm text-gray-600 uppercase tracking-widest">
                Molecule Viewer
              </p>
              <p className="text-xs text-gray-1500 mt-2">
                Enter a SMILES string to begin
              </p>
            </div>
          </div>
        )}

        {!loading && !error && result && (
          <div className="space-y-6">
            {/* Source Molecule and Patterns Side by Side */}
            <div className="flex flex-col lg:flex-row gap-6 items-start">
              {/* SVG Visualization - Larger box on left */}
              <div className="flex-shrink-0 w-full lg:w-auto">
                <h3 className="text-sm font-light mb-3 uppercase tracking-wider text-black">
                  Source Molecule
                </h3>
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
                    <p className="text-[9px] uppercase tracking-wider text-gray-600 mb-1.5">Queries</p>
                    <div className="space-y-0.5 mb-3">
                      {buildExpressionLines().map((line, idx) => (
                        <div key={idx} className="flex items-baseline gap-1.5">
                          <span className="text-sm font-medium text-gray-1000 w-6 text-left flex-shrink-0">{idx + 1}.</span>
                          <span className="text-sm font-mono font-medium text-black flex-1">{line.expression}</span>
                          <span className="text-[10px] text-gray-600 font-mono flex-shrink-0">({line.variants})</span>
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
                    <div className="flex items-center justify-between mb-4">
                      <h3 className="text-sm font-light uppercase tracking-wider text-black">
                        Patterns <span className="text-gray-600">({filteredMatches.length})</span>
                      </h3>
                      <div className="flex items-center gap-2.5">
                        {/* Variant count input */}
                        <div className="flex items-center gap-2 bg-gray-50 border border-gray-200 rounded px-2.5 py-1.5">
                          <label htmlFor="variant-count" className="text-[10px] uppercase tracking-wider text-gray-600 font-medium">
                            Variants
                          </label>
                          <input
                            id="variant-count"
                            type="number"
                            min="1"
                            max="100"
                            value={variantCount}
                            onChange={(e) => setVariantCount(Math.max(1, Math.min(100, parseInt(e.target.value) || 1)))}
                            className="w-12 bg-white border border-gray-300 rounded px-1.5 py-0.5 text-[11px] font-mono text-center focus:outline-none focus:border-black transition-colors"
                          />
                        </div>
                        
                        <button
                          onClick={handleAddQuery}
                          disabled={selectedPatterns.size === 0}
                          className={`w-7 h-7 rounded-full font-semibold text-base transition-all flex items-center justify-center ${
                            selectedPatterns.size > 0
                              ? 'bg-black text-white hover:bg-gray-800'
                              : 'bg-gray-200 text-gray-400 cursor-not-allowed'
                          }`}
                          title="Add pattern to queries"
                        >
                          +
                        </button>
                        {expressionGroups.length > 0 && (
                          <button
                            onClick={handleClearExpression}
                            className="w-7 h-7 rounded-full bg-gray-200 hover:bg-gray-300 text-sm text-gray-600 flex items-center justify-center transition-colors"
                            title="Clear queries"
                          >
                            ✕
                          </button>
                        )}
                      </div>
                    </div>
                    
                    {/* Pattern grid */}
                    <div className="grid gap-2" style={{ gridTemplateColumns: 'repeat(auto-fill, minmax(120px, 1fr))' }}>
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
                          <p className="text-[10px] font-medium mb-1 pl-6 leading-tight line-clamp-2 min-h-[24px]">{match.trivial_name.name}</p>
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

            {/* Bio-isostere Results */}
            {bioisostereResults.size > 0 && result && (() => {
              const filteredMatches = result.matches.filter(match => 
                match.trivial_name.group.toLowerCase().includes('cyclic') ||
                match.trivial_name.group.toLowerCase().includes('biological') ||
                match.trivial_name.group.toLowerCase().includes('ring')
              );

              return (
                <div className="mt-12">
                  <h2 className="text-lg font-light mb-6 uppercase tracking-wider text-black">
                    Bio-isostere Results
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
                                <span className="text-[10px] font-medium text-black leading-tight line-clamp-3 overflow-hidden">{originalPattern.trivial_name.name}</span>
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
                              
                              return (
                                <div
                                  key={matchIdx}
                                  className="border border-gray-300 rounded p-2 bg-white hover:border-gray-400 transition-colors flex flex-col w-full"
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
                                </div>
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
