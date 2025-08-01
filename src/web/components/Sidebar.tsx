import React from "react";

interface SidebarProps {
  isOpen: boolean;
  onClose: () => void;
}

export const Sidebar: React.FC<SidebarProps> = ({ isOpen, onClose }) => {
  return (
    <>
      {/* Backdrop/Overlay */}
      {isOpen && (
        <div 
          className="sidebar-backdrop"
          onClick={onClose}
          aria-label="Close sidebar"
        />
      )}
      
      {/* Sidebar Panel */}
      <aside className={`sidebar ${isOpen ? 'open' : ''}`}>
        <div className="sidebar-content">
          {/* Empty content placeholder as requested */}
          <div className="sidebar-placeholder">
            <p>Sidebar content will be added here</p>
          </div>
        </div>
      </aside>
    </>
  );
};