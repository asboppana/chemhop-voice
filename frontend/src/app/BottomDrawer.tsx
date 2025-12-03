import React, { useState, useEffect, useRef } from 'react';
import { motion, AnimatePresence, useMotionValue, useTransform, animate } from 'framer-motion';
import { useDrawer, type DrawerState } from '@/contexts/DrawerContext';
import { useChatContext } from '@/contexts/ChatContext';
import { useIsMobile } from '@/hooks/use-mobile';
import { ChatPanel } from '@/features/chat/ChatPanel';

// Configuration for drawer positioning and behavior
const DRAWER_CONFIG = {
  EXPANDED_Y: 70,
  HALF_Y_PERCENT: 0.50,
  COLLAPSED_VISIBLE_HEIGHT: 100,
  SCROLL_THRESHOLD: 8,         // Ignore tiny gestures
  REQUIRED_UPWARD_SCROLLS: 5,  // Require more distinct upward scroll events to prevent accidental collapse
};

// Visual spacing constants (in px)
const GUTTER_PX = 4; // Left/right outer margin for the drawer
const GAP_PX = 4; // Visible gap between drawer and chat (4px from each side)
const CHAT_WIDTH_PERCENT = 0.2; // Keep existing behavior (w-1/5)
const ANIM_MS = 350; // Shared animation duration for drawer+chat

interface BottomDrawerProps {
  routeKey: string | undefined;
  disablePulling?: boolean;
}

export const BottomDrawer: React.FC<BottomDrawerProps> = ({ routeKey, disablePulling: propDisablePulling = false }) => {
  const isMobile = useIsMobile();
  const { drawerState, setDrawerState, isDrawerVisible, drawerCards, disablePulling: ctxDisablePulling } = useDrawer();
  const pullingDisabled = propDisablePulling || ctxDisablePulling;
  const { 
    shouldShowInlineChat, 
    isDesktop, 
    setMode,
    state: { isOpen, mode }
  } = useChatContext();
  
  const [windowHeight, setWindowHeight] = useState(window.innerHeight);
  const [windowWidth, setWindowWidth] = useState(window.innerWidth);
  const [isDragging, setIsDragging] = useState(false);
  const dragStartY = useRef<number>(0);
  const dragStartState = useRef<DrawerState>('half');
  const contentRef = useRef<HTMLDivElement | null>(null);
  const atTopRef = React.useRef(true);
  const isAutoScrollingRef = React.useRef(false);
  const scrollElRef = React.useRef<HTMLElement | null>(null);
  // intent detection refs
  const gestureStartsAtTopRef = React.useRef(false);
  const lastWheelTsRef = React.useRef(0);
  const upwardScrollCountRef = React.useRef(0);
  const wasAtTopRef = React.useRef(true); // Track if we were at top in previous scroll event
  const reachedTopTimestampRef = React.useRef<number>(0); // Timestamp when we first reached the top
  const COLLAPSE_COOLDOWN_MS = 500; // Require 0.5s pause after reaching top before allowing collapse
  
  useEffect(() => {
    // Track both height and width to size & position drawer/chat precisely
    const handleResize = () => {
      setWindowHeight(window.innerHeight);
      setWindowWidth(window.innerWidth);
    };
    window.addEventListener('resize', handleResize);
    return () => window.removeEventListener('resize', handleResize);
  }, []);
  
  // Find primary scrollable element within drawer and track whether it's at the top
  function findPrimaryScrollElement(root: HTMLElement): HTMLElement {
    if (root.scrollHeight > root.clientHeight + 1) return root;
    let best: HTMLElement | null = null;
    try {
      const walker = document.createTreeWalker(root, NodeFilter.SHOW_ELEMENT);
      while (walker.nextNode()) {
        const node = walker.currentNode as HTMLElement;
        const style = window.getComputedStyle(node);
        const canScroll = style.overflowY !== 'hidden' && node.scrollHeight > node.clientHeight + 1;
        if (canScroll) {
          if (!best || node.scrollHeight > best.scrollHeight) best = node;
        }
      }
    } catch { /* no-op */ }
    return best || root;
  }

  useEffect(() => {
    const root = contentRef.current;
    if (!root) return;
    
    let cleanupScroll: (() => void) | null = null;
    let resizeObserver: ResizeObserver | null = null;
    let mutationObserver: MutationObserver | null = null;
    let retryTimeoutId: number | null = null;
    let retries = 0;

    const attach = () => {
      const primary = findPrimaryScrollElement(root);
      scrollElRef.current = primary;
      
      const onScroll = () => {
        const scrollTop = Math.max(0, primary.scrollTop);
        const isAtTop = scrollTop <= 1; // Allow 1px tolerance

        wasAtTopRef.current = atTopRef.current; // Store previous state
        
        // Track when we first reach the top
        if (!atTopRef.current && isAtTop) {
          // Just reached the top - record timestamp
          reachedTopTimestampRef.current = performance.now();
        } else if (!isAtTop) {
          // Scrolled away from top - reset timestamp
          reachedTopTimestampRef.current = 0;
        }
        
        atTopRef.current = isAtTop;
      };

      // Cleanup any previous listener before attaching a new one
      if (cleanupScroll) cleanupScroll();
      primary.addEventListener('scroll', onScroll, { passive: true });
      cleanupScroll = () => primary.removeEventListener('scroll', onScroll);

      // Initialize top state
      const initialScrollTop = Math.max(0, primary.scrollTop);
      const initialAtTop = initialScrollTop <= 1;
      atTopRef.current = initialAtTop;
      wasAtTopRef.current = initialAtTop;
      // If starting at top, set timestamp to allow immediate collapse attempts (already been at top)
      if (initialAtTop) {
        reachedTopTimestampRef.current = performance.now() - COLLAPSE_COOLDOWN_MS;
      }
    };

    const tryAttach = () => {
      attach();
      const primary = scrollElRef.current || root;
      // If we only found the root and it's not scrollable yet, retry briefly
      if (primary === root && !(primary.scrollHeight > primary.clientHeight + 1) && retries < 8) {
        retries += 1;
        retryTimeoutId = window.setTimeout(tryAttach, 50);
      }
    };

    // Initial pass after a frame so children render
    const rafId = window.requestAnimationFrame(() => { tryAttach(); });

    // Recalculate when size or DOM changes (content often mounts after first render)
    if ('ResizeObserver' in window) {
      resizeObserver = new ResizeObserver(() => { tryAttach(); });
      resizeObserver.observe(root);
    }
    mutationObserver = new MutationObserver(() => { tryAttach(); });
    mutationObserver.observe(root, { childList: true, subtree: true });

    return () => {
      if (cleanupScroll) cleanupScroll();
      if (resizeObserver) resizeObserver.disconnect();
      if (mutationObserver) mutationObserver.disconnect();
      if (retryTimeoutId) window.clearTimeout(retryTimeoutId);
      try { window.cancelAnimationFrame(rafId); } catch { /* no-op */ }
    };
  }, [routeKey]);

  function getWheelDeltaY(e: Pick<WheelEvent | React.WheelEvent, 'deltaY' | 'deltaMode'>): number {
    const dy = e.deltaMode === 1
      ? e.deltaY * 16
      : (e.deltaMode === 2 ? e.deltaY * window.innerHeight : e.deltaY);
    // Clamp extreme spikes from certain mice/browsers
    return Math.max(Math.min(dy, 400), -400);
  }

  function getAllScrollableElements(root: HTMLElement): HTMLElement[] {
    const result: HTMLElement[] = [];
    const walker = document.createTreeWalker(root, NodeFilter.SHOW_ELEMENT);
    while (walker.nextNode()) {
      const node = walker.currentNode as HTMLElement;
      const style = window.getComputedStyle(node);
      const canScroll = style.overflowY !== 'hidden' && node.scrollHeight > node.clientHeight + 1;
      if (canScroll) result.push(node);
    }
    // Include root last in case it is scrollable
    if (root.scrollHeight > root.clientHeight + 1) result.push(root);
    return result;
  }

  function areAllAtTop(elements: HTMLElement[]): boolean {
    for (const el of elements) {
      const scrollTop = Math.max(0, el.scrollTop); // Fix negative values
      if (scrollTop > 1) return false; // Allow 1px tolerance
    }
    return true;
  }

  function waitForAllAtTop(elements: HTMLElement[], maxFrames: number = 8): Promise<void> {
    return new Promise((resolve) => {
      let frames = 0;
      const check = () => {
        frames += 1;
        if (areAllAtTop(elements) || frames >= maxFrames) {
          return resolve();
        }
        requestAnimationFrame(check);
      };
      requestAnimationFrame(check);
    });
  }

  async function requestCollapseToHalf() {
    // Prevent multiple simultaneous collapse requests
    if (isAutoScrollingRef.current) return;
    
    const root = scrollElRef.current || contentRef.current;
    if (root) {
      isAutoScrollingRef.current = true;
      const scrollables = getAllScrollableElements(root);
      
      for (const s of scrollables) {
        const scrollTop = Math.max(0, s.scrollTop);
        if (scrollTop > 1) { // Use tolerance
          s.scrollTo({ top: 0, behavior: 'smooth' }); // Use smooth for better UX
        }
      }
      
      await waitForAllAtTop(scrollables, 15); // Increase timeout for smooth scroll
      atTopRef.current = true;
      isAutoScrollingRef.current = false;
    }
    
    setDrawerState('half');
  }
  
  const drawerPositions = {
    expanded: DRAWER_CONFIG.EXPANDED_Y,
    half: windowHeight * DRAWER_CONFIG.HALF_Y_PERCENT,
    collapsed: windowHeight - DRAWER_CONFIG.COLLAPSED_VISIBLE_HEIGHT,
  };
  // Use a shared motion value for vertical position so inner content height tracks the drawer smoothly
  const drawerY = useMotionValue<number>(drawerPositions[drawerState]);

  useEffect(() => {
    // Sync motion to target position; no animation when dragging
    const controls = animate(
      drawerY,
      drawerPositions[drawerState],
      { duration: isDragging ? 0 : ANIM_MS / 1000, ease: 'linear' }
    );
    return () => { try { controls.stop(); } catch { void 0; } };
  }, [drawerState, windowHeight, isDragging, drawerPositions, drawerY]);

  // Derived heights from the animated y-value
  const contentHeight = useTransform(drawerY, (y) => Math.max(0, windowHeight - y - 32));
  const chatHeight = useTransform(drawerY, (y) => Math.max(0, windowHeight - y));

  // Whether inline chat is actively visible on desktop
  const chatOpen = shouldShowInlineChat && isDesktop;

  // Pixel width of chat (kept proportional to viewport width)
  const chatWidthPx = chatOpen ? Math.round(windowWidth * CHAT_WIDTH_PERCENT) : 0;
  
  useEffect(() => {
    // Allow chatting with the drawer in half state; only bump from collapsed -> half
    if (!pullingDisabled && isOpen && mode === 'attached' && !isMobile) {
      if (drawerState === 'collapsed') setDrawerState('half');
    }
  }, [pullingDisabled, isOpen, mode, isMobile, setDrawerState, drawerState]);
  
  // Force detached mode when collapsed, but allow free choice in other states
  useEffect(() => {
    if (!pullingDisabled && isOpen && !isMobile) {
      if (drawerState === 'collapsed' && mode === 'attached') {
        setMode('detached');
      }
    }
  }, [pullingDisabled, drawerState, isOpen, isMobile, setMode, mode]);

  // Force expanded when pulling is disabled
  useEffect(() => {
    if (pullingDisabled && drawerState !== 'expanded') {
      setDrawerState('expanded');
    }
  }, [pullingDisabled, drawerState, setDrawerState]);

const handleDrawerBarClick = async (e: React.MouseEvent) => {
  if (pullingDisabled) return;
  e.stopPropagation();
  
  const getNextState = (current: DrawerState): DrawerState => {
    switch (current) {
      case 'collapsed': return 'half';
      case 'half': return 'expanded';
      case 'expanded': return 'half';
      default: return 'half';
    }
  };
  const next = getNextState(drawerState);

  // Ensure content is at top before collapsing from expanded -> half
  if (drawerState === 'expanded' && next === 'half') {
    const primaryElement = scrollElRef.current;
    if (primaryElement) {
      const scrollTop = Math.max(0, primaryElement.scrollTop);
      const isAtTop = scrollTop <= 1;
      if (!isAtTop) {
        await requestCollapseToHalf();
        return; // prevent collapse until content is at the top
      }
    }
  }

  setDrawerState(next);
};
  const handleTouchStart = (e: React.TouchEvent) => {
    if (pullingDisabled) return;
    setIsDragging(true);
    dragStartY.current = e.touches[0].clientY;
    dragStartState.current = drawerState;
  };

  const handleMouseDown = (e: React.MouseEvent) => {
    if (pullingDisabled) return;
    setIsDragging(true);
    dragStartY.current = e.clientY;
    dragStartState.current = drawerState;
  };

  const handleTouchMove = (e: React.TouchEvent) => {
    if (!isDragging || isAutoScrollingRef.current) return;

    const currentY = e.touches[0].clientY;
    const deltaY = currentY - dragStartY.current;
    const threshold = 50; // Minimum drag distance to trigger state change
    
    if (Math.abs(deltaY) > threshold) {
      const dragDirection = deltaY > 0 ? 'down' : 'up';
      let newState: DrawerState = dragStartState.current;
      
      if (dragDirection === 'up') {
        // Dragging up - expand
        if (dragStartState.current === 'collapsed') newState = 'half';
        else if (dragStartState.current === 'half') newState = 'expanded';
      } else {
        // Dragging down - collapse
        if (dragStartState.current === 'expanded') {
          const primaryElement = scrollElRef.current;
          if (primaryElement) {
            const scrollTop = Math.max(0, primaryElement.scrollTop);
            const isAtTop = scrollTop <= 1;
            if (!isAtTop) {
              return;
            }
          }
          newState = 'half';
        }
        else if (dragStartState.current === 'half') newState = 'collapsed';
      }
      
      if (newState !== drawerState) {
        setDrawerState(newState);
        setIsDragging(false); // End drag after state change
      }
    }
  };

  const handleMouseMove = (e: MouseEvent) => {
    if (!isDragging || isAutoScrollingRef.current) return;

    const currentY = e.clientY;
    const deltaY = currentY - dragStartY.current;
    const threshold = 50;
    
    if (Math.abs(deltaY) > threshold) {
      const dragDirection = deltaY > 0 ? 'down' : 'up';
      let newState: DrawerState = dragStartState.current;
      
      if (dragDirection === 'up') {
        if (dragStartState.current === 'collapsed') newState = 'half';
        else if (dragStartState.current === 'half') newState = 'expanded';
      } else {
        if (dragStartState.current === 'expanded') {
          const primaryElement = scrollElRef.current;
          if (primaryElement) {
            const scrollTop = Math.max(0, primaryElement.scrollTop);
            const isAtTop = scrollTop <= 1;
            if (!isAtTop) {
              return;
            }
          }
          newState = 'half';
        }
        else if (dragStartState.current === 'half') newState = 'collapsed';
      }
      
      if (newState !== drawerState) {
        setDrawerState(newState);
        setIsDragging(false);
      }
    }
  };

  const handleDragEnd = () => {
    setIsDragging(false);
  };

  // Add mouse move and up listeners when dragging
  useEffect(() => {
    if (!pullingDisabled && isDragging) {
      document.addEventListener('mousemove', handleMouseMove);
      document.addEventListener('mouseup', handleDragEnd);
      document.addEventListener('touchend', handleDragEnd);
      
      return () => {
        document.removeEventListener('mousemove', handleMouseMove);
        document.removeEventListener('mouseup', handleDragEnd);
        document.removeEventListener('touchend', handleDragEnd);
      };
    }
  }, [pullingDisabled, isDragging]);
  
  // Expand when scrolling down at top; collapse when scrolling up at top (desktop/trackpad)
  const handleContentWheel = React.useCallback(
    (e: React.WheelEvent<HTMLDivElement>) => {
      if (pullingDisabled || isDragging || isAutoScrollingRef.current) return;
      if (Math.abs(e.deltaX) > Math.abs(e.deltaY)) return; // ignore horizontal gestures
      
      // Ignore wheel events within interactive overlay components (dropdowns, popovers, dialogs)
      const target = e.target as HTMLElement;
      const isInOverlay = target.closest('[role="dialog"], [role="menu"], [role="listbox"], [data-radix-popper-content-wrapper]');
      if (isInOverlay) return;

      const now = e.timeStamp || performance.now();
      const dt = now - lastWheelTsRef.current;
      lastWheelTsRef.current = now;

      // --- ROBUST "AT TOP" CHECK ---
      // Instead of relying on a single cached 'scrollElRef', dynamically walk up 
      // from the event target to check if ANY scrollable ancestor is scrolled down.
      let isDeep = false;
      let current: HTMLElement | null = target;
      const root = contentRef.current;
      
      // Safety break to prevent infinite loops in detached subtrees
      let depth = 0;
      while (current && root && root.contains(current) && depth < 50) {
        // Optimization: Only check computed style if scroll height indicates it might be scrollable
        if (current.scrollHeight > current.clientHeight + 1) {
           const style = window.getComputedStyle(current);
           if (style.overflowY === 'auto' || style.overflowY === 'scroll') {
               // Found a scrollable container. Is it scrolled down?
               if (current.scrollTop > 1) {
                   isDeep = true;
                   break;
               }
           }
        }
        if (current === root) break;
        current = current.parentElement;
        depth++;
      }

      const isCurrentlyAtTop = !isDeep;
      
      // Store previous top state from this wheel handler's perspective
      const wasPreviouslyAtTop = wasAtTopRef.current;

      // If significant time passed, it's a new gesture - reset everything
      if (dt > 200) {
        gestureStartsAtTopRef.current = isCurrentlyAtTop;
        upwardScrollCountRef.current = 0;
        wasAtTopRef.current = isCurrentlyAtTop;
        
        // If we are starting a new gesture at the top, ensure we have a timestamp
        // But if we were ALREADY at the top (persisted), keep the old timestamp to allow immediate action
        if (isCurrentlyAtTop) {
            if (reachedTopTimestampRef.current === 0) {
                reachedTopTimestampRef.current = now;
            }
        } else {
            reachedTopTimestampRef.current = 0;
        }
      }

      const dy = getWheelDeltaY(e);
      const { SCROLL_THRESHOLD, REQUIRED_UPWARD_SCROLLS } = DRAWER_CONFIG;

      const scrollingDown = dy > SCROLL_THRESHOLD;
      const scrollingUp = dy < -SCROLL_THRESHOLD;

      // Reset counter on any downward scroll
      if (scrollingDown) {
        upwardScrollCountRef.current = 0;
        // If scrolling down away from top, reset the timestamp
        if (!isCurrentlyAtTop) {
          reachedTopTimestampRef.current = 0;
        }
      }

      // From collapsed/half -> expanded on downward scroll
      if (scrollingDown && isCurrentlyAtTop && drawerState !== 'expanded') {
        setDrawerState('expanded');
        return;
      }

      // From expanded -> half on upward scroll: intent-based with counter
      if (scrollingUp && drawerState === 'expanded') {
        // 1. If we are NOT at the top (we are deep), we definitely shouldn't collapse.
        if (!isCurrentlyAtTop) {
          upwardScrollCountRef.current = 0;
          gestureStartsAtTopRef.current = false;
          wasAtTopRef.current = false;
          reachedTopTimestampRef.current = 0;
          return;
        }

        // 2. If we just reached the top in this specific event...
        if (!wasPreviouslyAtTop && isCurrentlyAtTop) {
          reachedTopTimestampRef.current = now;
          upwardScrollCountRef.current = 0;
          gestureStartsAtTopRef.current = true;
          wasAtTopRef.current = true;
          return;
        }

        // 3. Enforce cooldown: must wait at least 0.5s after reaching top before collapse
        const timeSinceReachedTop = now - reachedTopTimestampRef.current;
        if (reachedTopTimestampRef.current > 0 && timeSinceReachedTop < COLLAPSE_COOLDOWN_MS) {
          upwardScrollCountRef.current = 0;
          wasAtTopRef.current = true;
          return;
        }

        // 4. Must have started this gesture at the top
        if (!gestureStartsAtTopRef.current) {
          upwardScrollCountRef.current = 0;
          wasAtTopRef.current = true;
          return;
        }

        // 5. Finally, we can count this scroll toward collapse
        upwardScrollCountRef.current += 1;
        wasAtTopRef.current = true;
        if (upwardScrollCountRef.current >= REQUIRED_UPWARD_SCROLLS) {
          upwardScrollCountRef.current = 0;
          setDrawerState('half');
        }
      }
    },
    [pullingDisabled, isDragging, drawerState, setDrawerState]
  );

  if (!isDrawerVisible) return null;
  
  // We switch from width classes to left/right offsets to:
  // - add 4px gutters on both sides
  // - animate the drawer being "pushed" left by the chat via the right offset
  // Right offset is 4px when chat is closed, and (chat width + 8px gap) when open
  const drawerRightOffset = chatOpen ? chatWidthPx + GAP_PX : GUTTER_PX;
  
  return (
    <>
      <motion.div
        className={`fixed bottom-0 bg-white z-30 rounded-t-[24px] overflow-visible`}
        style={{ 
          left: GUTTER_PX,
          height: windowHeight,
          y: drawerY,
        }}
        // Animate horizontal push (right) and vertical position (y) together
        initial={{ right: drawerRightOffset }}
        animate={{ right: drawerRightOffset }}
        transition={{ duration: isDragging ? 0 : ANIM_MS / 1000, ease: 'linear' }}
        onTouchStart={handleTouchStart}
        onTouchMove={handleTouchMove}
        onMouseDown={handleMouseDown}
      >
        <div
          className={`w-full h-8 flex items-center justify-center ${pullingDisabled ? 'cursor-default' : 'cursor-pointer'} rounded-t-[24px]`}
          onClick={pullingDisabled ? undefined : handleDrawerBarClick}
        >
          {!pullingDisabled && (
            <div className="w-16 h-1 bg-gray-500 hover:bg-gray-1000 rounded-full" />
          )}
        </div>

        <motion.div 
          className="bg-white overflow-y-auto"
          ref={contentRef}
          onWheel={handleContentWheel}
          style={{ height: contentHeight }}
        >
          <AnimatePresence mode="wait">
            <motion.div
              key={routeKey}
              initial={{ opacity: 0, y: 20 }}
              animate={{ opacity: 1, y: 0 }}
              exit={{ opacity: 0, y: -20 }}
              transition={{ duration: 0.2, ease: 'linear' }}
              className="h-full"
            >
              <div className="flex flex-row gap-6 w-full h-full px-6 pt-2 pb-0">
                {drawerCards.map((card, index) => (
                  <div key={index} className={`${card.width} min-w-[260px] h-full`}>
                    {card.content}
                  </div>
                ))}
              </div>
            </motion.div>
          </AnimatePresence>
        </motion.div>
      </motion.div>

      <AnimatePresence>
        {chatOpen && (
          <motion.div
            className={`fixed bottom-0 right-0 z-30`}
            // Slide in horizontally; height tracks the drawer so top aligns exactly with drawer's top
            initial={{ x: '100%' }}
            animate={{ x: 0 }}
            exit={{ x: '100%' }}
            transition={{ duration: isDragging ? 0 : ANIM_MS / 1000, ease: 'linear' }}
            style={{ width: chatWidthPx, height: chatHeight }}
          >
            <ChatPanel isInlineMode={true} parentControlsInlineAnimation={true} />
          </motion.div>
        )}
      </AnimatePresence>
    </>
  );
};