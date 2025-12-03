import React, { createContext, useContext, useState, useEffect, type ReactNode, useCallback, useMemo, useRef } from 'react';
import { useLocation } from 'react-router-dom';

import { type DrawerCard } from '@/app/PageWithDrawerLayout';
import { getAllRoutes } from '@/app/routes';

export type DrawerState = 'collapsed' | 'half' | 'expanded';

interface DrawerContextType {
  drawerState: DrawerState;
  setDrawerState: (state: DrawerState) => void;
  resetToHalf: () => void;
  drawerCards: DrawerCard[];
  setDrawerCards: (cards: DrawerCard[]) => void;
  isDrawerVisible: boolean;
  setIsDrawerVisible: (visible: boolean) => void;
  initialDrawerState: DrawerState;
  setInitialDrawerState: (state: DrawerState) => void;
  setManualDrawerState: (state: DrawerState | null) => void;
  disablePulling: boolean;
  setDisablePulling: (disabled: boolean) => void;
}

// Separate context for actions only (to prevent unnecessary re-renders)
interface DrawerActionsContextType {
  setDrawerState: (state: DrawerState) => void;
  resetToHalf: () => void;
  setDrawerCards: (cards: DrawerCard[]) => void;
  setIsDrawerVisible: (visible: boolean) => void;
  setInitialDrawerState: (state: DrawerState) => void;
  setManualDrawerState: (state: DrawerState | null) => void;
  setDisablePulling: (disabled: boolean) => void;
}

const DrawerContext = createContext<DrawerContextType | null>(null);
const DrawerActionsContext = createContext<DrawerActionsContextType | null>(null);

/**
 * Determines the appropriate drawer state based on the current URL path
 * 
 * This function implements URL-based drawer positioning rules:
 * - Specific detail pages (e.g., /results/systems/cardiovascular) → expanded
 * - General section pages (e.g., /results/systems) → half
 * - Other pages → use fallback state
 */
export const getDrawerStateFromPath = (pathname: string, fallbackState: DrawerState = 'half'): DrawerState => {
  // 1) Path-specific rules for /results take precedence over route defaults
  //    - /results and /results/ → half
  //    - /results/systems/summary → half
  //    - Any other /results/* → expanded
  // if (pathname === '/results' || pathname === '/results/') {
  //   return 'half';
  // }
  if (pathname.match(/^\/(results\/systems\/summary|results\/systems\/?$)\/?$/)) {
    return 'half';
  }
  // if (pathname.startsWith('/results/')) {
  //   return 'expanded';
  // }

  // 2) Otherwise, look up the route configuration to get the proper initialDrawerState
  const allRoutes = getAllRoutes();
  const currentRoute = allRoutes.find(route => {
    // Check for exact match
    if (pathname === route.path) {
      return true;
    }
    
    // Check for wildcard matches
    if (route.path.includes('*')) {
      const basePath = route.path.replace('/*', '');
      return pathname.startsWith(basePath);
    }
    
    // Check for parametrized routes
    if (route.path.includes(':')) {
      const routeParts = route.path.split('/');
      const locationParts = pathname.split('/');
      
      if (routeParts.length !== locationParts.length) {
        return false;
      }
      
      return routeParts.every((part, index) => 
        part.startsWith(':') || part === locationParts[index]
      );
    }
    
    return false;
  });

  // If route disables pulling, force expanded regardless of any other settings
  if (currentRoute?.layoutConfig?.disableDrawerPulling) {
    return 'expanded';
  }
  if (currentRoute?.layoutConfig?.initialDrawerState) {
    return currentRoute.layoutConfig.initialDrawerState;
  }

  // 3) Fallback legacy rules
  if (pathname.match(/^\/results\/systems\/[^\/]+$/) && pathname !== '/results/systems/summary') {
    return 'expanded';
  }
  if (pathname.match(/^\/results\/records\/[^\/]+$/)) {
    return 'expanded';
  }
  if (pathname.match(/^\/results\/biomarkers\/[^\/]+$/)) {
    return 'expanded';
  }
  if (pathname.startsWith('/results')) {
    return 'half';
  }

  // 4) Final fallback
  return fallbackState;
};

/**
 * Compute a stable "route key" for a pathname based on the matching route definition.
 * For wildcard routes like "/results/*", we return the base path ("/results").
 * Returns null if no matching route definition is found.
 */
const getRouteKeyForPath = (pathname: string): string | null => {
  const allRoutes = getAllRoutes();
  const currentRoute = allRoutes.find(route => {
    // Exact match
    if (pathname === route.path) return true;

    // Wildcard match
    if (route.path.includes('*')) {
      const basePath = route.path.replace('/*', '');
      return pathname.startsWith(basePath);
    }

    // Parametrized match
    if (route.path.includes(':')) {
      const routeParts = route.path.split('/');
      const locationParts = pathname.split('/');
      if (routeParts.length !== locationParts.length) return false;
      return routeParts.every((part, index) => part.startsWith(':') || part === locationParts[index]);
    }

    return false;
  });

  if (!currentRoute) return null;
  return currentRoute.path.includes('/*') ? currentRoute.path.split('/*')[0] : currentRoute.path;
};

/**
 * DrawerProvider - Context provider for managing global drawer state
 * 
 * This provider implements a sophisticated state management system:
 * 1. URL-based default states (automatic positioning based on route)
 * 2. Manual overrides (component-triggered state changes)
 * 3. Automatic cleanup (manual overrides cleared on navigation)
 * 4. Optimized re-rendering (actions context separate from state context)
 * 
 * This ensures consistent behavior while allowing components to control the drawer when needed.
 */
export const DrawerProvider: React.FC<{children: ReactNode}> = ({ children }) => {
  const location = useLocation();
  
  const [isDrawerVisible, setIsDrawerVisible] = useState(false);
  const [drawerCards, setDrawerCards] = useState<DrawerCard[]>([]);
  const [initialDrawerState, setInitialDrawerState] = useState<DrawerState>('half');
  const [disablePulling, setDisablePulling] = useState<boolean>(false);
  
  // Track the previous pathname to detect navigation
  const [previousPathname, setPreviousPathname] = useState(location.pathname);
  // Track the previous route key (parent/base route) to reset overrides across sections
  const [previousRouteKey, setPreviousRouteKey] = useState<string | null>(getRouteKeyForPath(location.pathname));
  
  // Manual state override - allows components to temporarily control drawer position
  const [manualDrawerState, setManualDrawerState] = useState<DrawerState | null>(null);

  // Calculate the default state based on the current URL
  const urlBasedState = getDrawerStateFromPath(location.pathname, initialDrawerState);
  
  // The final, effective state prioritizes manual overrides over URL-based defaults
  const drawerState = manualDrawerState ?? urlBasedState;

  // Only clear manual overrides when the URL-based state would change
  useEffect(() => {
    // Only proceed if pathname actually changed
    if (location.pathname !== previousPathname) {
      const previousUrlBasedState = getDrawerStateFromPath(previousPathname, initialDrawerState);
      const newUrlBasedState = getDrawerStateFromPath(location.pathname, initialDrawerState);

      const newRouteKey = getRouteKeyForPath(location.pathname);

      // Clear manual override when either:
      // 1) URL-based state would change, or
      // 2) The base route key (parent route) changed (e.g., /dashboard -> /results)
      if (previousUrlBasedState !== newUrlBasedState || previousRouteKey !== newRouteKey) {
        setManualDrawerState(null);
      }

      setPreviousPathname(location.pathname);
      setPreviousRouteKey(newRouteKey);
    }
  }, [location.pathname, previousPathname, previousRouteKey, initialDrawerState]);

  // Memoize the functions to prevent unnecessary re-renders
  const setDrawerStateAndOverride = useCallback((state: DrawerState) => {
    setManualDrawerState(state);
  }, []);
  
  const resetToHalf = useCallback(() => setManualDrawerState('half'), []);

  // Use ref to store actions so they don't cause re-renders
  const actionsRef = useRef<DrawerActionsContextType>({
    setDrawerState: setDrawerStateAndOverride,
    resetToHalf,
    setDrawerCards,
    setIsDrawerVisible,
    setInitialDrawerState,
    setManualDrawerState,
    setDisablePulling
  });

  // Update the ref when functions change (but this won't cause re-renders)
  actionsRef.current = {
    setDrawerState: setDrawerStateAndOverride,
    resetToHalf,
    setDrawerCards,
    setIsDrawerVisible,
    setInitialDrawerState,
    setManualDrawerState,
    setDisablePulling
  };

  // Actions context value - stable reference, never changes
  const actionsContextValue = useMemo(() => actionsRef.current, []);

  // Main context value - includes state that will cause re-renders
  const contextValue = useMemo(() => ({ 
    drawerState, 
    setDrawerState: setDrawerStateAndOverride, 
    resetToHalf,
    drawerCards, 
    setDrawerCards,
    isDrawerVisible, 
    setIsDrawerVisible,
    initialDrawerState, 
    setInitialDrawerState,
    setManualDrawerState,
    disablePulling,
    setDisablePulling
  }), [
    drawerState, 
    setDrawerStateAndOverride, 
    resetToHalf, 
    drawerCards, 
    isDrawerVisible, 
    initialDrawerState,
    disablePulling
  ]);
  
  return (
    <DrawerContext.Provider value={contextValue}>
      <DrawerActionsContext.Provider value={actionsContextValue}>
        {children}
      </DrawerActionsContext.Provider>
    </DrawerContext.Provider>
  );
};

/**
 * Hook for accessing drawer context (includes state - will cause re-renders)
 * Use this for components that need to read drawer state
 */
export const useDrawer = () => {
  const context = useContext(DrawerContext);
  if (!context) {
    console.warn('useDrawer used outside DrawerProvider - returning fallback state');
    // Return fallback state instead of throwing
    return {
      drawerState: 'half' as DrawerState,
      setDrawerState: () => {},
      resetToHalf: () => {},
      drawerCards: [],
      setDrawerCards: () => {},
      isDrawerVisible: false,
      setIsDrawerVisible: () => {},
      initialDrawerState: 'half' as DrawerState,
      setInitialDrawerState: () => {},
      setManualDrawerState: () => {},
      disablePulling: false,
      setDisablePulling: () => {}
    };
  }
  return context;
};

/**
 * Hook for accessing only drawer actions (no state - won't cause re-renders)
 * Use this for components that only need to control the drawer, not read its state
 */
export const useDrawerActions = () => {
  const actionsContext = useContext(DrawerActionsContext);
  if (actionsContext) {
    return actionsContext;
  }

  // Fallback: if actions context is missing but main context exists, derive actions from it
  const stateContext = useContext(DrawerContext);
  if (stateContext) {
    const {
      setDrawerState,
      resetToHalf,
      setDrawerCards,
      setIsDrawerVisible,
      setInitialDrawerState,
      setManualDrawerState,
      setDisablePulling,
    } = stateContext;
    return {
      setDrawerState,
      resetToHalf,
      setDrawerCards,
      setIsDrawerVisible,
      setInitialDrawerState,
      setManualDrawerState,
      setDisablePulling,
    } as DrawerActionsContextType;
  }

  console.warn('useDrawerActions used outside DrawerProvider - returning no-op functions');
  // Return no-op functions instead of throwing
  return {
    setDrawerState: () => {},
    resetToHalf: () => {},
    setDrawerCards: () => {},
    setIsDrawerVisible: () => {},
    setInitialDrawerState: () => {},
    setManualDrawerState: () => {},
    setDisablePulling: () => {}
  };
}; 