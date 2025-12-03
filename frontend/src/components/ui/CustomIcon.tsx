import React from 'react';

interface CustomIconProps {
  name: string;
  size?: number;
  className?: string;
  color?: string;
}

export const CustomIcon: React.FC<CustomIconProps> = ({ 
  name, 
  size = 20, 
  className = "", 
  color = "currentColor" 
}) => {
  return (
    <img 
      src={`/icons/${name}.svg`}
      alt={name}
      width={size}
      height={size}
      className={className}
      style={{ 
        filter: color === "currentColor" 
          ? 'brightness(0) saturate(100%) invert(100%)' // Makes black SVGs white
          : `brightness(0) saturate(100%) ${color}`
      }}
    />
  );
};

// Fixed approach for dynamic imports
export const CustomIconOptimized: React.FC<CustomIconProps> = ({ 
  name, 
  size = 20, 
  className = "" 
}) => {
  const [iconSrc, setIconSrc] = React.useState<string>('');
  const [loading, setLoading] = React.useState<boolean>(true);

  React.useEffect(() => {
    const loadIcon = async () => {
      try {
        // Use the public URL approach instead of dynamic import
        const iconUrl = `/icons/${name}.svg`;
        
        // Verify the icon exists
        const response = await fetch(iconUrl);
        if (response.ok) {
          setIconSrc(iconUrl);
        } else {
          console.warn(`Icon ${name} not found`);
        }
      } catch (error) {
        console.warn(`Error loading icon ${name}:`, error);
      } finally {
        setLoading(false);
      }
    };
    
    loadIcon();
  }, [name]);

  if (loading) {
    return (
      <div 
        style={{ width: size, height: size }} 
        className={className}
      />
    );
  }

  if (!iconSrc) return null;

  return (
    <img 
      src={iconSrc}
      alt={name}
      width={size}
      height={size}
      className={className}
    />
  );
};

// Alternative approach using a more Vite-friendly pattern
export const CustomIconDynamic: React.FC<CustomIconProps> = ({ 
  name, 
  size = 20, 
  className = "" 
}) => {
  const [iconContent, setIconContent] = React.useState<string>('');

  React.useEffect(() => {
    const loadIcon = async () => {
      try {
        // This approach works better with Vite
        const iconUrl = new URL(`../../../public/icons/${name}.svg`, import.meta.url).href;
        const response = await fetch(iconUrl);
        if (response.ok) {
          setIconContent(await response.text());
        }
      } catch (error) {
        console.warn(`Error loading icon ${name}:`, error);
        // Fallback to simple URL approach
        setIconContent(`<img src="/icons/${name}.svg" alt="${name}" />`);
      }
    };
    
    loadIcon();
  }, [name]);

  if (!iconContent) return null;

  return (
    <div 
      className={className}
      style={{ width: size, height: size }}
      dangerouslySetInnerHTML={{ __html: iconContent }}
    />
  );
}; 