import React, { useState, useRef, useEffect } from 'react';
import { Link, useLocation } from 'react-router-dom';
import { motion, AnimatePresence } from 'framer-motion';

import { CustomIcon } from '@/components/ui/CustomIcon';
import HoverRotatingSedonaLogo from '@/components/animations/HoverRotatingSedonaLogo';

interface HeaderProps {
  userRole?: string;
}

const Header: React.FC<HeaderProps> = () => {
  const location = useLocation();
  const [isProfileDropdownOpen, setProfileDropdownOpen] = useState(false);
  const dropdownRef = useRef<HTMLDivElement>(null);
  const homeRef = useRef<HTMLAnchorElement | null>(null);

  const [indicatorLeft, setIndicatorLeft] = useState(0);
  const [indicatorWidth, setIndicatorWidth] = useState(0);

  // Close dropdown when clicking outside
  useEffect(() => {
    const handleClickOutside = (event: MouseEvent) => {
      if (dropdownRef.current && !dropdownRef.current.contains(event.target as Node)) {
        setProfileDropdownOpen(false);
      }
    };

    document.addEventListener("mousedown", handleClickOutside);
    return () => {
      document.removeEventListener("mousedown", handleClickOutside);
    };
  }, []);

  // Update navigation indicator position for Home link
  useEffect(() => {
    const isHomeActive = location.pathname === '/' || location.pathname === '/home';
    if (isHomeActive && homeRef.current) {
      const { offsetLeft, offsetWidth } = homeRef.current;
      setIndicatorLeft(offsetLeft);
      setIndicatorWidth(offsetWidth);
    } else {
      setIndicatorWidth(0);
    }
  }, [location.pathname]);

  const isHomeActive = location.pathname === '/' || location.pathname === '/home';

  const ProfileCircle = ({ textSize = 'text-xs', isHoverable = false }) => {
    const size = 'w-8 h-8';
    const initials = 'U'; // Default initials for public user
    
    const circleContent = (
      <span className={`${textSize} font-medium text-white`}>
        {initials}
      </span>
    );

    if (isHoverable) {
      return (
        <motion.div
          className={`${size} rounded-full bg-gray-1500 flex items-center justify-center overflow-hidden`}
          whileHover={{ scale: 1.08 }}
          transition={{ duration: 0.3 }}
        >
          {circleContent}
        </motion.div>
      );
    }

    return (
      <div className={`${size} bg-gray-1500 rounded-full flex items-center justify-center overflow-hidden`}>
        {circleContent}
      </div>
    );
  };

  return (
    <motion.div 
      initial={{ opacity: 0, y: -20 }}
      animate={{ opacity: 1, y: 0 }}
      transition={{ duration: 0.6 }}
      className="fixed top-0 left-0 right-0 z-[60] bg-transparent mb-4"
    >
      <div className="grid grid-cols-[1fr,auto,1fr] items-center px-6 pt-2 pb-2">
        {/* Left Side - Logo */}
        <div className="flex items-center space-x-4 justify-self-start">
          <HoverRotatingSedonaLogo />
        </div>

        {/* Center - Home Navigation */}
        <div className="relative flex items-center space-x-4 justify-self-center">
          <div 
            className="absolute bg-white/10 rounded-full transition-all duration-300 ease-in-out"
            style={{ 
              left: `${indicatorLeft}px`, 
              width: `${indicatorWidth}px`,
              height: '36px',
              opacity: indicatorWidth > 0 ? 1 : 0
            }}
          />
          <Link
            ref={homeRef}
            to="/"
            className={`
              group flex uppercase items-center space-x-2 px-3 py-2 rounded-full transition-all duration-200
              !no-underline focus:!outline-none 
              ${isHomeActive 
                ? '!text-white' 
                : '!text-white/70 hover:!text-white hover:bg-white/10 hover:shadow-gray-1000'
              }
              visited:!text-white focus:!text-white
            `}
            title="Home"
          >
            <CustomIcon 
              name="Home" 
              size={18} 
              className="text-current transition-transform duration-200"
            />
            <span className="text-[13px] font-medium">Home</span>
          </Link>
        </div>

        {/* Right Side - Profile Dropdown */}
        <div className="flex items-center space-x-2 justify-self-end">
          <div
            className="relative"
            ref={dropdownRef}
            onMouseEnter={() => setProfileDropdownOpen(true)}
            onMouseLeave={() => setProfileDropdownOpen(false)}
          >
            <motion.button
              className="flex items-center space-x-2 p-2 text-white/70 hover:text-white rounded-xl focus:outline-none"
              title="Profile"
              whileHover="hover"
              initial="rest"
              tabIndex={0}
              aria-haspopup="true"
              aria-expanded={isProfileDropdownOpen}
            >
              <ProfileCircle textSize="text-xs" isHoverable={true} />
            </motion.button>
            
            <AnimatePresence>
              {isProfileDropdownOpen && (
                <motion.div
                  initial={{ opacity: 0, scale: 0.95, y: -10, backgroundColor: 'rgba(0,0,0,0.35)' }}
                  animate={{ 
                    opacity: 1, 
                    scale: 1, 
                    y: 0,
                    backgroundColor: isProfileDropdownOpen ? 'rgba(0,0,0,0.8)' : 'rgba(0,0,0,0.35)'
                  }}
                  exit={{ opacity: 0, scale: 0.95, y: -10, backgroundColor: 'rgba(0,0,0,0.35)' }}
                  transition={{ duration: 0.15 }}
                  className="absolute right-0 mt-0 w-48 origin-top-right rounded-md backdrop-blur-lg shadow-lg"
                  style={{ zIndex: 100 }}
                >
                  <div className="py-0" role="menu" aria-orientation="vertical" aria-labelledby="user-menu-button">
                    <button 
                      onClick={() => {
                        setProfileDropdownOpen(false);
                        // Add logout logic here
                        console.log('Logout clicked');
                      }} 
                      className="flex items-center px-4 py-2 text-sm text-white/80 hover:bg-white/10 hover:text-white hover:rounded-md w-full text-left" 
                      role="menuitem"
                    >
                      <CustomIcon name="Logout" size={16} className="mr-2" />
                      <span>Logout</span>
                    </button>
                  </div>
                </motion.div>
              )}
            </AnimatePresence>
          </div>
        </div>
      </div>
    </motion.div>
  );
};

export default Header;