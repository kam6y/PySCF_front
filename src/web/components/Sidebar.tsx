import React from "react";
import { CalculationInstance } from "../types/calculation";

interface SidebarProps {
  isOpen: boolean;
  onClose: () => void;
  calculations: CalculationInstance[];
  activeCalculationId: string | null;
  calculationsLoading: boolean;
  calculationsError: string | null;
  onCalculationSelect: (calculationId: string) => void;
  onCalculationRename: (calculationId: string, newName: string) => Promise<void>;
  onCalculationDelete: (calculationId: string) => Promise<void>;
}

export const Sidebar: React.FC<SidebarProps> = ({ 
  isOpen, 
  onClose, 
  calculations, 
  activeCalculationId, 
  calculationsLoading, 
  calculationsError, 
  onCalculationSelect, 
  onCalculationRename, 
  onCalculationDelete 
}) => {
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
          <div className="sidebar-header">
            <h3>Calculations</h3>
          </div>
          
          <div className="sidebar-calculations">
            {calculationsLoading ? (
              <div className="sidebar-loading">
                <p>Loading calculations...</p>
              </div>
            ) : calculationsError ? (
              <div className="sidebar-error">
                <p>Error: {calculationsError}</p>
              </div>
            ) : calculations.length === 0 ? (
              <div className="sidebar-empty">
                <p>No calculations yet</p>
                <p>Click the + button to create one</p>
              </div>
            ) : (
              calculations.map((calculation) => (
                <div 
                  key={calculation.id}
                  className={`sidebar-calculation-item ${
                    calculation.id === activeCalculationId ? 'active' : ''
                  }`}
                  onClick={() => onCalculationSelect(calculation.id)}
                >
                  <div className="calculation-info">
                    <div className="calculation-name">{calculation.name}</div>
                    <div className="calculation-status">
                      <span className={`status-badge ${calculation.status}`}>
                        {calculation.status}
                      </span>
                    </div>
                  </div>
                  <div className="calculation-actions">
                    <button 
                      className="rename-btn"
                      onClick={(e) => {
                        e.stopPropagation();
                        const newName = prompt('Enter new name:', calculation.name);
                        if (newName && newName !== calculation.name) {
                          onCalculationRename(calculation.id, newName);
                        }
                      }}
                    >
                      ✏️
                    </button>
                    <button 
                      className="delete-btn"
                      onClick={(e) => {
                        e.stopPropagation();
                        if (confirm(`Delete calculation "${calculation.name}"?`)) {
                          onCalculationDelete(calculation.id);
                        }
                      }}
                    >
                      🗑️
                    </button>
                  </div>
                </div>
              ))
            )}
          </div>
        </div>
      </aside>
    </>
  );
};