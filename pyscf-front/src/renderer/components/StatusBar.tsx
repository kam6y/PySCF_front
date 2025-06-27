import React, { useState } from 'react'

interface StatusBarProps {
  isCalculating: boolean
  logs: Array<{ type: string; data: string }>
  onClearLogs: () => void
}

const StatusBar: React.FC<StatusBarProps> = ({ 
  isCalculating, 
  logs, 
  onClearLogs 
}) => {
  const [showLogs, setShowLogs] = useState(false)

  const getStatusColor = () => {
    if (isCalculating) return 'text-warning'
    if (logs.some(log => log.type === 'error' || log.type === 'stderr')) {
      return 'text-error'
    }
    return 'text-success'
  }

  const getStatusText = () => {
    if (isCalculating) return 'Calculating...'
    if (logs.some(log => log.type === 'error' || log.type === 'stderr')) {
      return 'Error'
    }
    return 'Ready'
  }

  return (
    <>
      <div className="h-8 bg-tertiary border-t border-border flex items-center px-4 text-sm">
        <div className="flex items-center gap-2">
          <div className={`w-2 h-2 rounded-full ${
            isCalculating ? 'bg-warning animate-pulse' : 
            logs.some(log => log.type === 'error') ? 'bg-error' : 'bg-success'
          }`} />
          <span className={getStatusColor()}>
            {getStatusText()}
          </span>
        </div>

        <div className="flex-1" />

        <div className="flex items-center gap-4 text-xs text-secondary">
          <span>PySCF Backend Connected</span>
          <button
            onClick={() => setShowLogs(!showLogs)}
            className="hover:text-primary"
          >
            {logs.length > 0 && (
              <span className="mr-1">({logs.length})</span>
            )}
            Logs {showLogs ? '▼' : '▲'}
          </button>
        </div>
      </div>

      {/* Log Panel */}
      {showLogs && (
        <div className="h-32 bg-primary border-t border-border flex flex-col">
          <div className="flex items-center justify-between px-4 py-2 bg-secondary border-b border-border">
            <span className="text-sm font-medium">System Logs</span>
            <div className="flex gap-2">
              <button
                onClick={onClearLogs}
                className="text-xs text-secondary hover:text-primary"
              >
                Clear
              </button>
              <button
                onClick={() => setShowLogs(false)}
                className="text-xs text-secondary hover:text-primary"
              >
                ✕
              </button>
            </div>
          </div>
          
          <div className="flex-1 overflow-y-auto p-2 font-mono text-xs">
            {logs.length === 0 ? (
              <div className="text-secondary italic">No logs yet...</div>
            ) : (
              logs.map((log, index) => (
                <div
                  key={index}
                  className={`mb-1 ${
                    log.type === 'error' || log.type === 'stderr'
                      ? 'text-error'
                      : log.type === 'stdout'
                      ? 'text-info'
                      : 'text-secondary'
                  }`}
                >
                  <span className="text-muted">
                    [{new Date().toLocaleTimeString()}]
                  </span>{' '}
                  {log.data}
                </div>
              ))
            )}
          </div>
        </div>
      )}
    </>
  )
}

export default StatusBar