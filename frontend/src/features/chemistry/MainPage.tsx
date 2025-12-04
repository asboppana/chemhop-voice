import React, { useState } from 'react';
import { analyzeMolecule } from '@/services/moleculeService';
import type { MoleculeAnalysisResponse } from '@/services/moleculeService';
import { SmartsRenderer } from '@/components/SmartsRenderer';

export const ChemistryMainPage: React.FC = () => {
  const [smilesInput, setSmilesInput] = useState('');
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [result, setResult] = useState<MoleculeAnalysisResponse | null>(null);

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
      console.log('Molecule analysis result:', response);
    } catch (err: any) {
      setError(err.response?.data?.detail || 'Failed to analyze molecule');
      console.error('Error analyzing molecule:', err);
    } finally {
      setLoading(false);
    }
  };

  return (
    <div className="h-screen w-full bg-black text-white flex">
      {/* Left side - SMILES Input */}
      <div className="w-80 border-r border-gray-800 flex flex-col p-6">
        <h1 className="text-xl font-light mb-6 tracking-wide">SMILES INPUT</h1>
        
        <form onSubmit={handleSubmit} className="flex flex-col gap-4">
          <div className="flex flex-col gap-2">
            <label htmlFor="smiles" className="text-xs uppercase tracking-wider text-gray-400">
              Enter SMILES String
            </label>
            <input
              id="smiles"
              type="text"
              value={smilesInput}
              onChange={(e) => setSmilesInput(e.target.value)}
              placeholder="CCO"
              className="bg-transparent border border-gray-700 rounded px-3 py-2 text-sm focus:outline-none focus:border-white transition-colors font-mono"
            />
          </div>
          
          <button
            type="submit"
            className="bg-white text-black px-4 py-2 text-sm font-medium hover:bg-gray-200 transition-colors"
          >
            SUBMIT
          </button>
        </form>

        {/* Example SMILES */}
        <div className="mt-8 pt-6 border-t border-gray-800">
          <p className="text-xs uppercase tracking-wider text-gray-400 mb-3">Examples</p>
          <div className="flex flex-col gap-2 max-h-96 overflow-y-auto">
            {[
              { name: 'Sildenafil (Viagra)', smiles: 'CCCc1nn(C)c2c1nc(nc2=O)c3c(OCC)cccc3S(=O)(=O)N4CCN(C)CC4' },
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
                className="text-left p-2 hover:bg-gray-900 transition-colors rounded text-xs"
              >
                <div className="font-medium">{example.name}</div>
                <div className="text-gray-500 font-mono text-[10px] mt-1">{example.smiles}</div>
              </button>
            ))}
          </div>
        </div>
      </div>

      {/* Right side - Main content area */}
      <div className="flex-1 p-12 overflow-auto">
        {loading && (
          <div className="flex items-center justify-center h-full">
            <div className="text-center">
              <div className="animate-spin rounded-full h-12 w-12 border-b-2 border-white mx-auto mb-4"></div>
              <p className="text-sm text-gray-400">Analyzing molecule...</p>
            </div>
          </div>
        )}

        {error && (
          <div className="flex items-center justify-center h-full">
            <div className="text-center max-w-md">
              <p className="text-sm text-red-400 mb-2">Error</p>
              <p className="text-xs text-gray-400">{error}</p>
            </div>
          </div>
        )}

        {!loading && !error && !result && (
          <div className="flex items-center justify-center h-full">
            <div className="text-center">
              <p className="text-sm text-gray-500 uppercase tracking-widest">
                Molecule Viewer
              </p>
              <p className="text-xs text-gray-600 mt-2">
                Enter a SMILES string to begin
              </p>
            </div>
          </div>
        )}

        {!loading && !error && result && (
          <div className="space-y-4">
            {/* Source Molecule and Patterns Side by Side */}
            <div className="flex gap-6 items-start">
              {/* SVG Visualization - Larger box on left */}
              <div className="flex-shrink-0">
                <h3 className="text-sm font-light mb-3 uppercase tracking-wider text-gray-400">
                  Source Molecule
                </h3>
                <div className="border border-gray-800 rounded-lg p-4 bg-white overflow-hidden" style={{ width: '280px', height: '280px' }}>
                  <div 
                    dangerouslySetInnerHTML={{ __html: result.svg }}
                    className="w-full h-full"
                    style={{ 
                      display: 'flex',
                      alignItems: 'center',
                      justifyContent: 'center'
                    }}
                  />
                </div>
                <p className="text-[10px] text-gray-600 font-mono mt-2">{result.smiles}</p>
              </div>

              {/* Matched Patterns as mini boxes - filter for cyclic/biological */}
              {result.matches && result.matches.length > 0 && (() => {
                const filteredMatches = result.matches.filter(match => 
                  match.trivial_name.group.toLowerCase().includes('cyclic') ||
                  match.trivial_name.group.toLowerCase().includes('biological') ||
                  match.trivial_name.group.toLowerCase().includes('ring')
                );
                
                return filteredMatches.length > 0 ? (
                  <div className="flex-1">
                    <h3 className="text-sm font-light mb-3 uppercase tracking-wider text-gray-400">
                      Patterns ({filteredMatches.length})
                    </h3>
                    <div className="grid grid-cols-3 gap-3 auto-rows-max">
                      {filteredMatches.map((match, index) => (
                        <div 
                          key={index} 
                          className="border border-gray-800 rounded p-2 hover:border-gray-700 transition-colors bg-gray-900 flex flex-col"
                          style={{ width: '180px', height: '180px' }}
                        >
                          <p className="text-xs font-medium mb-1 truncate">{match.trivial_name.name}</p>
                          <p className="text-[9px] text-gray-500 mb-2 truncate">{match.trivial_name.group}</p>
                          <div className="flex-1 overflow-hidden">
                            <SmartsRenderer smarts={match.trivial_name.smarts} />
                          </div>
                        </div>
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
