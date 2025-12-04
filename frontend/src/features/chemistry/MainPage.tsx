import React, { useState } from 'react';
import { analyzeMolecule, highlightMolecule } from '@/services/moleculeService';
import type { MoleculeAnalysisResponse } from '@/services/moleculeService';

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
  const [expressionGroups, setExpressionGroups] = useState<number[][]>([]); // Array of AND groups, ORed together

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

  // Toggle pattern selection (multi-select)
  const handlePatternClick = async (_atomIndices: number[], patternIndex: number) => {
    if (!result) return;
    
    const newSelected = new Set(selectedPatterns);
    if (newSelected.has(patternIndex)) {
      newSelected.delete(patternIndex);
    } else {
      newSelected.add(patternIndex);
    }
    setSelectedPatterns(newSelected);
    
    // Highlight all selected patterns
    if (newSelected.size > 0) {
      // Get filtered matches to access atom indices
      const filteredMatches = result.matches.filter(match => 
        match.trivial_name.group.toLowerCase().includes('cyclic') ||
        match.trivial_name.group.toLowerCase().includes('biological') ||
        match.trivial_name.group.toLowerCase().includes('ring')
      );
      
      // Combine atom indices from all selected patterns
      const allAtomIndices: number[] = [];
      newSelected.forEach(idx => {
        if (filteredMatches[idx]) {
          allAtomIndices.push(...filteredMatches[idx].atom_indices);
        }
      });
      
      // Remove duplicates
      const uniqueAtomIndices = [...new Set(allAtomIndices)];
      
      setHighlightLoading(true);
      try {
        const response = await highlightMolecule(analyzedSmiles, uniqueAtomIndices);
        setDisplayedSvg(response.svg);
      } catch (err: any) {
        console.error('Error highlighting molecule:', err);
      } finally {
        setHighlightLoading(false);
      }
    } else {
      setDisplayedSvg(result.svg);
    }
  };

  // Handle AND button - group selected patterns together
  const handleAndClick = () => {
    if (selectedPatterns.size === 0) return;
    
    const selectedArray = Array.from(selectedPatterns).sort((a, b) => a - b);
    
    // Add to existing expression
    if (expressionGroups.length === 0) {
      // First group
      setExpressionGroups([selectedArray]);
    } else {
      // Add to the last OR group (extend with AND)
      const newGroups = [...expressionGroups];
      const lastGroup = newGroups[newGroups.length - 1];
      // Merge with last group (AND operation within same OR group)
      const mergedGroup = [...new Set([...lastGroup, ...selectedArray])].sort((a, b) => a - b);
      newGroups[newGroups.length - 1] = mergedGroup;
      setExpressionGroups(newGroups);
    }
    
    // Clear selection
    setSelectedPatterns(new Set());
  };

  // Handle OR button - start a new group
  const handleOrClick = () => {
    if (selectedPatterns.size === 0) return;
    
    const selectedArray = Array.from(selectedPatterns).sort((a, b) => a - b);
    
    // Add as a new OR group
    setExpressionGroups([...expressionGroups, selectedArray]);
    
    // Clear selection
    setSelectedPatterns(new Set());
  };

  // Build expression groups for display (each OR group on its own line)
  const buildExpressionLines = (): string[] => {
    if (expressionGroups.length === 0) return [];
    
    return expressionGroups.map(group => 
      group.map(idx => indexToLetter(idx)).join(' AND ')
    );
  };

  // Clear expression
  const handleClearExpression = () => {
    setExpressionGroups([]);
    setSelectedPatterns(new Set());
    if (result) {
      setDisplayedSvg(result.svg);
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
                <div className="border border-gray-300 rounded-lg p-4 bg-white w-full max-w-sm mx-auto lg:mx-0 aspect-square lg:w-80 relative">
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
                <p className="text-[10px] text-gray-600 font-mono mt-2 break-all max-w-sm mx-auto lg:mx-0">{result.smiles}</p>
                
                {/* Expression display */}
                {expressionGroups.length > 0 && (
                  <div className="mt-4 p-3 bg-gray-100 rounded-lg max-w-sm mx-auto lg:mx-0">
                    <p className="text-[9px] uppercase tracking-wider text-gray-600 mb-1.5">Expression</p>
                    <div className="space-y-0.5">
                      {buildExpressionLines().map((line, idx) => (
                        <div key={idx} className="flex items-baseline gap-1.5">
                          <span className="text-sm font-medium text-gray-1000 w-6 text-left flex-shrink-0">{idx + 1}.</span>
                          <span className="text-sm font-mono font-medium text-black">{line}</span>
                          {idx < buildExpressionLines().length - 1 && (
                            <span className="text-sm font-semibold text-gray-1000">OR</span>
                          )}
                        </div>
                      ))}
                    </div>
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
                    {/* Header with AND/OR buttons */}
                    <div className="flex items-center justify-between mb-3">
                      <h3 className="text-sm font-light uppercase tracking-wider text-black">
                        Patterns ({filteredMatches.length})
                      </h3>
                      <div className="flex items-center gap-2">
                        <button
                          onClick={handleAndClick}
                          disabled={selectedPatterns.size === 0}
                          className={`px-3 py-1 text-[11px] font-semibold rounded border-2 transition-all ${
                            selectedPatterns.size > 0
                              ? 'border-black bg-black text-white hover:bg-gray-800 hover:border-gray-800 shadow-sm'
                              : 'border-gray-300 bg-gray-100 text-gray-400 cursor-not-allowed'
                          }`}
                        >
                          AND
                        </button>
                        <button
                          onClick={handleOrClick}
                          disabled={selectedPatterns.size === 0}
                          className={`px-3 py-1 text-[11px] font-semibold rounded border-2 transition-all ${
                            selectedPatterns.size > 0
                              ? 'border-black bg-white text-black hover:bg-gray-50 shadow-sm'
                              : 'border-gray-300 bg-gray-100 text-gray-400 cursor-not-allowed'
                          }`}
                        >
                          OR
                        </button>
                        {expressionGroups.length > 0 && (
                          <button
                            onClick={handleClearExpression}
                            className="w-6 h-6 rounded-full bg-gray-300 text-xs text-gray-600 flex items-center justify-center"
                            title="Clear expression"
                          >
                            ✕
                          </button>
                        )}
                      </div>
                    </div>
                    
                    {/* Pattern grid */}
                    <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-3 gap-3">
                      {filteredMatches.map((match, index) => (
                        <button 
                          key={index}
                          onClick={() => handlePatternClick(match.atom_indices, index)}
                          className={`border rounded p-2 transition-colors flex flex-col aspect-square cursor-pointer text-left relative max-w-[160px] ${
                            selectedPatterns.has(index) 
                              ? 'border-blue-500 bg-blue-50 ring-2 ring-blue-200' 
                              : 'border-gray-300 bg-gray-50 hover:border-gray-400'
                          }`}
                        >
                          {/* Letter label badge */}
                          <div className="absolute top-1 left-1 w-5 h-5 rounded-full bg-black text-white text-[10px] font-medium flex items-center justify-center">
                            {indexToLetter(index)}
                          </div>
                          <p className="text-xs font-medium mb-2 truncate pl-6">{match.trivial_name.name}</p>
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
          </div>
        )}
      </div>
    </div>
  );
};

export default ChemistryMainPage;
