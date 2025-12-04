import React, { useEffect, useState } from 'react';
import { analyzeMolecule } from '@/services/moleculeService';

interface SmartsRendererProps {
  smarts: string;
}

export const SmartsRenderer: React.FC<SmartsRendererProps> = ({ smarts }) => {
  const [svg, setSvg] = useState<string | null>(null);
  const [loading, setLoading] = useState(true);

  useEffect(() => {
    const fetchSvg = async () => {
      try {
        // Try to analyze the SMARTS pattern as if it were SMILES
        const result = await analyzeMolecule(smarts);
        setSvg(result.svg);
      } catch (error) {
        // If it fails, just show the SMARTS text
        console.error('Failed to render SMARTS:', error);
      } finally {
        setLoading(false);
      }
    };

    fetchSvg();
  }, [smarts]);

  if (loading) {
    return (
      <div className="bg-white rounded h-full flex items-center justify-center">
        <div className="w-4 h-4 border-2 border-gray-300 border-t-gray-600 rounded-full animate-spin"></div>
      </div>
    );
  }

  if (!svg) {
    return (
      <div className="bg-white rounded h-full flex items-center justify-center p-2">
        <p className="text-[8px] text-gray-400 text-center break-all leading-tight">
          {smarts}
        </p>
      </div>
    );
  }

  return (
    <div className="bg-white rounded h-full w-full overflow-hidden flex items-center justify-center">
      <div 
        dangerouslySetInnerHTML={{ __html: svg }}
        style={{ 
          maxWidth: '100%',
          maxHeight: '100%',
          display: 'flex',
          alignItems: 'center',
          justifyContent: 'center'
        }}
      />
    </div>
  );
};
