import React from 'react';

/**
 * SIMPLIFIED ROUTING SYSTEM
 * 
 * Three core routes using PageWithDrawerLayout:
 * - /home: Main dashboard
 * - /chat: Chat interface with optional id parameter
 * - /profile: User profile
 */

import Placeholder from '@/app/Placeholder';
import ChemistryMainPage from '@/features/chemistry/MainPage';

// Core navigation configuration
export const MAIN_NAVIGATION = [
  { label: 'Home', key: 'HOME', path: '/home' },
  { label: 'Chat', key: 'CHAT', path: '/chat' },
  { label: 'Profile', key: 'PROFILE', path: '/profile' }
];

// Simplified route definition - all routes use page-with-drawer layout
export interface RouteDefinition {
  label: string;
  path: string;
  component: React.ComponentType<any>;
  isPublic?: boolean;
  layout: 'page-with-drawer';
  showChatBar?: boolean; // Optional: show chat bar (defaults to true)
  
  // Layout configuration for PageWithDrawerLayout
  layoutConfig?: {
    mainBackground?: 'default' | 'blur' | 'solid-black' | 'solid-white';
    initialDrawerState?: 'collapsed' | 'half' | 'expanded';
    drawerComponent?: React.ComponentType<any>;
    chatMode?: 'floating-only' | 'inline-allowed';
    disableDrawerPulling?: boolean;
  };
}

// Three core routes - all using page-with-drawer layout
export const routeDefinitions: RouteDefinition[] = [
  {
    label: "Home",
    path: "/home",
    component: Placeholder,
    layout: 'page-with-drawer',
    isPublic: true,
    layoutConfig: {
      mainBackground: 'default',
      initialDrawerState: 'expanded',
      drawerComponent: ChemistryMainPage,
      chatMode: 'inline-allowed',
      disableDrawerPulling: true
    }
  },
];

// ============================================================================
// HELPER FUNCTIONS
// ============================================================================

// Get all routes
export const getAllRoutes = (): RouteDefinition[] => {
  return routeDefinitions;
};


// Public routes that don't require authentication
export const publicRoutes = routeDefinitions.filter(route => route.isPublic).map(route => route.path);
