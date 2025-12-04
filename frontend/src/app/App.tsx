import React, { useState, useMemo } from "react";
import { BrowserRouter as Router, Route, Routes, Navigate, useLocation, Outlet } from "react-router-dom";
import { Toaster } from "@/components/ui/toaster";
import { useToast } from "@/hooks/use-toast";
import "./index.css";
import "./global.css";

import { Context } from "@/contexts/Context";
import { ChatProvider, useChatContext } from '@/contexts/ChatContext';
import { getAllRoutes, type RouteDefinition } from "./routes";
import { ChatPanel } from '@/features/chat/ChatPanel';
import { FloatingChatBubble } from '@/features/chat/FloatingChatBubble';
import { AskBar } from '@/features/chat/AskBar';
import Header from "@/app/Header";
import { PageWithDrawerLayout, type DrawerCard } from '@/app/PageWithDrawerLayout';
import { DrawerProvider } from '@/contexts/DrawerContext';
import { BottomDrawer } from '@/app/BottomDrawer';


// ============================================================================
// MAIN LAYOUT COMPONENT
// Orchestrates all persistent UI elements and animations
// ============================================================================
const MainLayout = () => {
	const location = useLocation();

	// Get chat context using new reducer pattern (moved early so effects can close chats)
	const chatCtx = useChatContext();
	const { 
		shouldShowInlineChat, 
		shouldShowFloatingChat,
		state: { mode, isOpen, messages },
	} = chatCtx;

	// Find current route definition to determine layout, background, etc.
	const allRoutes = getAllRoutes();
	const currentRoute = allRoutes.find(route => {
		// Check for exact match
		if (location.pathname === route.path) {
			return true;
		}
		
		// Check for wildcard matches
		if (route.path.includes('*')) {
			const basePath = route.path.replace('/*', '');
			return location.pathname.startsWith(basePath);
		}
		
		// Check for parametrized routes
		if (route.path.includes(':')) {
			const routeParts = route.path.split('/');
			const locationParts = location.pathname.split('/');
			
			// Handle routes with both parameters and wildcards
			// For routes like "/admin/users/:userId/results/*"
			if (route.path.includes('*')) {
				// Find the wildcard position
				const wildcardIndex = routeParts.findIndex(part => part === '*');
				if (wildcardIndex === -1) return false;
				
				// Check if we have at least as many location parts as route parts up to the wildcard
				if (locationParts.length < wildcardIndex) return false;
				
				// Check that all parts before the wildcard match (including parameters)
				for (let i = 0; i < wildcardIndex; i++) {
					const routePart = routeParts[i];
					const locationPart = locationParts[i];
					
					if (!routePart.startsWith(':') && routePart !== locationPart) {
						return false;
					}
				}
				
				return true;
			}
			
			// Handle regular parametrized routes (no wildcards)
			if (routeParts.length !== locationParts.length) {
				return false;
			}
			
			return routeParts.every((part, index) => 
				part.startsWith(':') || part === locationParts[index]
			);
		}
		
		return false;
	});

	const isPageWithDrawerLayout = currentRoute?.layout === 'page-with-drawer';
	const showUserHeader = isPageWithDrawerLayout;
	
	// Chat configuration - only for USER role
	const currentChatMode = currentRoute?.layoutConfig?.chatMode || 'inline-allowed';
	const showChatBar = currentRoute?.showChatBar !== false; // Default to true unless explicitly set to false
	const showAskBar = !isOpen && showChatBar;

	// Show floating bubble when chat is closed and conditions are met
	const hasMessages = messages.length > 0;
	const showFloatingBubble = !isOpen && (
		(showAskBar && mode === 'detached') || 
		hasMessages ||
		(showAskBar && currentChatMode === 'floating-only')
	);

	return (
		<>
			<div>
				{/* Navigation UI */}
				{showUserHeader && (
					<Header
						userRole={'USER'}
					/>
				)}

                {/* Main application shell */}
                <div className={`app-shell relative font-pp-neue overflow-hidden ${shouldShowInlineChat ? 'chat-open' : ''}`}>

					{/* Main content area with page transition animations */}
					<div className={`main-area relative transition-all duration-300 ${shouldShowInlineChat ? 'with-chat' : ''}`}>
                        <div className={`main-content`}>
							<Outlet />
						</div>
					</div>
				</div>
				
				{/* Persistent bottom drawer - only for page-with-drawer layout */}
				{isPageWithDrawerLayout && (
					<BottomDrawer 
						routeKey={currentRoute?.path} 
						disablePulling={currentRoute?.layoutConfig?.disableDrawerPulling === true}
					/>
				)}
			
				{/* Chat elements - only for USER role */}
				{shouldShowFloatingChat && (
					<ChatPanel isInlineMode={false} chatMode={currentChatMode} />
				)}

				{showFloatingBubble && <FloatingChatBubble />}
				{showAskBar && <AskBar />}
			</div>
		</>
	);
};

// ============================================================================
// NOT FOUND COMPONENT
// ============================================================================
const NotFound = () => {
	const { toast } = useToast();
	toast({
		title: "Page Not Found",
		description: "Redirecting you back home",
		variant: "destructive",
	});
	return <Navigate to={"/"} replace={true} />;
};

// ============================================================================
// LAYOUT WRAPPER COMPONENTS
// These components handle the different layout types
// ============================================================================

// Wrapper for page-with-drawer layout
const PageWithDrawerWrapper: React.FC<{ route: RouteDefinition }> = ({ route }) => {	
	// Use more stable memoization for desktop/drawer mode
	const DrawerComponent = route.layoutConfig?.drawerComponent;
	
	// Create a stable key based on the base route path
	// For /results/*, this will always be "/results" regardless of sub-path
	const baseRouteKey = route.path.split('/*')[0] || route.path;
	
	const drawerComponent = useMemo(() => 
		DrawerComponent ? <DrawerComponent /> : null,
		[DrawerComponent, baseRouteKey] // Stable for all sub-routes
	);
	
	const drawerCards: DrawerCard[] = useMemo(() => 
		drawerComponent ? [{
			content: drawerComponent,
			width: 'w-full' as const
		}] : [],
		[drawerComponent]
	);
	
	// Memoize the main content with base route key
	const mainContent = useMemo(() => 
		<route.component />,
		[route.component, baseRouteKey]
	);

	
	return (
		<PageWithDrawerLayout
			mainContent={mainContent}
			drawerCards={drawerCards}
			initialDrawerState={'half'}
		/>
	);
};


// ============================================================================
// ROOT APP COMPONENT
// Sets up the complete application architecture
// ============================================================================
const App: React.FC = () => {	
	const [context, setContext] = useState(null);
	const allRoutes = getAllRoutes();
	// Use flag to avoid double fading when overlay is active (checked inline where needed)
	
	// Route rendering logic - handles different layout types consistently
	const renderRoutes = (routes: RouteDefinition[]) => {
		return routes.map(route => {
			return (
			<Route 
				key={route.path} 
				path={route.path} 
				element={<PageWithDrawerWrapper route={route}/>} 
			/>
			);
		});
	};

	return (
		<Context.Provider value={[context, setContext]}>
				<ChatProvider>
					<Router>
						<DrawerProvider>
							<Routes>
								<Route path="/" element={<Navigate to="/home" replace />} />
								<Route element={
									<MainLayout />
								}>
									{renderRoutes(allRoutes)}
								</Route>
								<Route path="*" element={<NotFound />} />
							</Routes>
						</DrawerProvider>
						{/* Global black overlay for cross-route fade-in to avoid white flashes */}
						<Toaster />
					</Router>
				</ChatProvider>
		</Context.Provider>
	);
};

export default App;
