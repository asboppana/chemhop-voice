import React from 'react';
import { cn } from '@/lib/utils';

type XCloseSize = 'sm' | 'md';
type XCloseRenderAs = 'button' | 'icon';

interface XCloseProps {
  onClick?: React.MouseEventHandler;
  ariaLabel?: string;
  className?: string;
  size?: XCloseSize; // sm: 20px, md: 24px
  invertIcon?: boolean; // whether to invert the icon color
  renderAs?: XCloseRenderAs; // 'button' renders a button wrapper, 'icon' renders only the img
  borderColor?: string; // border color of the button
}

const sizeToButton: Record<XCloseSize, string> = {
  sm: 'w-5 h-5',
  md: 'w-6 h-6',
};

const sizeToIcon: Record<XCloseSize, string> = {
  sm: 'w-3.5 h-3.5',
  md: 'w-4 h-4',
};

export const XClose: React.FC<XCloseProps> = ({
  onClick,
  ariaLabel = 'Close',
  className,
  size = 'md',
  invertIcon = true,
  renderAs = 'button',
  borderColor = 'border-gray-1500',
}) => {
  if (renderAs === 'icon') {
    return (
      <img
        src="/icons/X.svg"
        alt={ariaLabel}
        className={cn(sizeToIcon[size], invertIcon && 'invert', className)}
        onClick={onClick}
      />
    );
  }
  return (
    <button
      type="button"
      onClick={onClick}
      aria-label={ariaLabel}
      className={cn(
        'rounded-full border flex items-center justify-center hover:bg-gray-1000 transition-colors',
        borderColor,
        sizeToButton[size],
        className,
      )}
    >
      <img
        src="/icons/X.svg"
        alt={ariaLabel}
        className={cn(sizeToIcon[size], invertIcon && 'invert')}
      />
    </button>
  );
};

export default XClose;


