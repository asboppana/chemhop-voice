import React, { useEffect, useCallback } from 'react';
import { useDrawerActions } from '@/contexts/DrawerContext';

export interface DrawerCard {
  content: React.ReactNode;
  width: 'w-3/5' | 'w-2/5' | 'w-1/5' | 'w-1/3' | 'w-2/3' | 'w-full';
}

interface MainContentLayoutProps {
  key?: string;
  mainContent?: React.ReactNode;
  drawerCards: DrawerCard[];
  initialDrawerState?: 'collapsed' | 'half' | 'expanded';
}

/**
 * PageWithDrawerLayout - Layout component for pages that include a bottom drawer
 * 
 * This component serves as the interface between individual pages and the global drawer system.
 * It handles:
 * - Configuring the persistent drawer with page-specific content
 * - Setting up the main content area above the drawer
 * - Managing pointer events when the drawer is expanded
 * 
 * This component is automatically wrapped around pages that specify 'page-with-drawer' layout
 * in their route configuration.
 */
export const PageWithDrawerLayout: React.FC<MainContentLayoutProps> = ({
  mainContent,
  drawerCards,
  initialDrawerState = 'half',
}) => {
  const { 
    setDrawerState,
    setDrawerCards, 
    setIsDrawerVisible, 
    setInitialDrawerState
  } = useDrawerActions();

  // Expand drawer when user scrolls up on page background; ignore mostly-horizontal gestures
  const handleBackgroundWheel = useCallback((e: React.WheelEvent<HTMLDivElement>) => {
    if (Math.abs(e.deltaX) > Math.abs(e.deltaY)) return;
    if (e.deltaY > 8) {
      setDrawerState('expanded');
    }
  }, [setDrawerState]);

  // Configure the persistent drawer when this page mounts
  useEffect(() => {
    // Show the drawer and set its content
    setIsDrawerVisible(true);
    setDrawerCards(drawerCards);
    setInitialDrawerState(initialDrawerState);

    // Clean up when page unmounts (hide drawer)
    return () => {
      setIsDrawerVisible(false);
      setDrawerCards([]);
    };
  }, [drawerCards, initialDrawerState, setDrawerCards, setIsDrawerVisible, setInitialDrawerState]);

  // Handle drawer card updates separately to avoid re-running the entire effect
  useEffect(() => {
    setDrawerCards(drawerCards);
  }, [drawerCards, setDrawerCards]);

  // Handle initial state updates separately
  useEffect(() => {
    setInitialDrawerState(initialDrawerState);
  }, [initialDrawerState, setInitialDrawerState]);

  return (
    <div className="h-auto md:h-[100dvh]">
      {/* Main Content Area - The space above the drawer where page content is rendered */}
      <div 
        className="relative z-[30] pt-20 pb-8 px-8 h-full overflow-y-auto md:overflow-hidden"
        onWheel={handleBackgroundWheel}
      >
        {mainContent}
      </div>
    </div>
  );
};

export default PageWithDrawerLayout;