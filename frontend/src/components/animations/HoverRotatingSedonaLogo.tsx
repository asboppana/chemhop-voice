import React from 'react';

interface HoverRotatingSedonaLogoProps {
  size?: number;
  className?: string;
  color?: "white" | "black";
}

const HoverRotatingSedonaLogo: React.FC<HoverRotatingSedonaLogoProps> = ({ 
  size = 30,
  className = "",
  color = "white"
}) => {
  // Heavy easing: cubic-bezier(0.8, 0, 0.2, 1) - Slow start, fast middle, slow end
  const transitionTimingFunction = 'cubic-bezier(0.8, 0, 0.2, 1)';

  return (
    <div className={`inline-flex items-center justify-center group ${className}`}>
      <svg 
        width={size} 
        height={size} 
        viewBox="0 0 56 57" 
        fill="none" 
        xmlns="http://www.w3.org/2000/svg"
      >
        {/* Bottom Left Quarter Circle */}
        <path 
          d="M13.9629 28.2686C21.6757 28.2686 27.9629 34.5558 27.9629 42.2686C27.9626 49.9811 21.6755 56.2686 13.9629 56.2686H13.5977V50.0537H13.9629C18.276 50.0537 21.7488 46.5451 21.749 42.2686C21.749 37.9918 18.2762 34.4824 13.9629 34.4824C9.64975 34.4826 6.21387 37.9919 6.21387 42.2686V42.6338H0V42.2686C0 34.5559 6.25027 28.2688 13.9629 28.2686Z" 
          fill={color === "white" ? "white" : "black"}
          className="transition-transform duration-700 group-hover:rotate-180"
          style={{
            transformOrigin: '13.96px 42.27px',
            transitionTimingFunction
          }}
        />
        
        {/* Bottom Right Quarter Circle */}
        <path 
          d="M41.9648 28.2686C49.6776 28.2686 55.9648 34.5558 55.9648 42.2686V42.6338H49.751V42.2686C49.751 37.9553 46.2416 34.4825 41.9648 34.4824C37.6881 34.4824 34.1787 37.9552 34.1787 42.2686C34.179 46.5817 37.6882 50.0176 41.9648 50.0176H42.3301V56.2314H41.9648C34.2522 56.2314 27.9651 49.9811 27.9648 42.2686C27.9648 34.5558 34.252 28.2686 41.9648 28.2686Z" 
          fill={color === "white" ? "white" : "black"}
          className="transition-transform duration-700 group-hover:rotate-180"
          style={{
            transformOrigin: '41.96px 42.27px',
            transitionTimingFunction
          }}
        />
        
        {/* Top Right Quarter Circle */}
        <path 
          d="M42.3672 6.48242H42.001C37.6878 6.48253 34.2148 9.99186 34.2148 14.2686C34.2151 18.545 37.6879 22.0536 42.001 22.0537C46.3141 22.0537 49.7507 18.5451 49.751 14.2686V13.9023H55.9648V14.2686C55.9646 21.9811 49.7136 28.2686 42.001 28.2686C34.2884 28.2684 28.0012 21.9811 28.001 14.2686C28.001 6.55583 34.2883 0.268666 42.001 0.268555H42.3672V6.48242Z" 
          fill={color === "white" ? "white" : "black"}
          className="transition-transform duration-700 group-hover:rotate-180"
          style={{
            transformOrigin: '42px 14.27px',
            transitionTimingFunction
          }}
        />
        
        {/* Top Left Quarter Circle */}
        <path 
          d="M14 0.286133C21.7128 0.286133 28 6.53721 28 14.25C27.9999 21.9627 21.7127 28.25 14 28.25C6.28728 28.25 0.000113415 21.9627 0 14.25V13.8838H6.21387V14.25C6.21398 18.5632 9.72331 22.0361 14 22.0361C18.2767 22.0361 21.786 18.5632 21.7861 14.25C21.7861 9.93668 18.2768 6.5 14 6.5H13.6348V0.286133H14Z" 
          fill={color === "white" ? "white" : "black"}
          className="transition-transform duration-700 group-hover:rotate-180"
          style={{
            transformOrigin: '14px 14.27px',
            transitionTimingFunction
          }}
        />
      </svg>
      

    </div>
  );
};

export default HoverRotatingSedonaLogo;