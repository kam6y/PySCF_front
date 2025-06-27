import React, { useState } from 'react'

interface CalculationPanelProps {
  onCalculate: (params: any) => void
  isCalculating: boolean
  results: any
}

const CalculationPanel: React.FC<CalculationPanelProps> = ({
  onCalculate,
  isCalculating,
  results
}) => {
  const [method, setMethod] = useState('DFT')
  const [functional, setFunctional] = useState('B3LYP')
  const [basis, setBasis] = useState('6-31G*')

  const handleSubmit = (e: React.FormEvent) => {
    e.preventDefault()
    onCalculate({
      method,
      functional: method === 'DFT' ? functional : undefined,
      basis
    })
  }

  const methods = ['HF', 'DFT', 'MP2']
  const functionals = ['B3LYP', 'PBE', 'PBE0', 'M06-2X', 'wB97X-D']
  const basisSets = ['STO-3G', '6-31G', '6-31G*', '6-31G**', '6-311G*', 'cc-pVDZ', 'cc-pVTZ']

  return (
    <div className="h-full flex flex-col">
      {/* Header */}
      <div className="p-4 border-b border-border">
        <h2 className="text-lg font-semibold">Quantum Chemistry Calculation</h2>
        <p className="text-sm text-secondary">Configure and run calculations</p>
      </div>

      <div className="flex-1 overflow-y-auto">
        {/* Calculation Parameters */}
        <div className="p-4">
          <form onSubmit={handleSubmit} className="space-y-4">
            {/* Method Selection */}
            <div>
              <label htmlFor="calculation-method" className="block text-sm font-medium mb-2">
                Calculation Method
              </label>
              <select
                id="calculation-method"
                value={method}
                onChange={(e) => setMethod(e.target.value)}
                className="input"
                disabled={isCalculating}
              >
                {methods.map(m => (
                  <option key={m} value={m}>{m}</option>
                ))}
              </select>
            </div>

            {/* Functional Selection (DFT only) */}
            {method === 'DFT' && (
              <div>
                <label htmlFor="exchange-functional" className="block text-sm font-medium mb-2">
                  Exchange-Correlation Functional
                </label>
                <select
                  id="exchange-functional"
                  value={functional}
                  onChange={(e) => setFunctional(e.target.value)}
                  className="input"
                  disabled={isCalculating}
                >
                  {functionals.map(f => (
                    <option key={f} value={f}>{f}</option>
                  ))}
                </select>
              </div>
            )}

            {/* Basis Set Selection */}
            <div>
              <label htmlFor="basis-set" className="block text-sm font-medium mb-2">
                Basis Set
              </label>
              <select
                id="basis-set"
                value={basis}
                onChange={(e) => setBasis(e.target.value)}
                className="input"
                disabled={isCalculating}
              >
                {basisSets.map(b => (
                  <option key={b} value={b}>{b}</option>
                ))}
              </select>
            </div>

            {/* Submit Button */}
            <button
              type="submit"
              disabled={isCalculating}
              className={`btn w-full ${
                isCalculating ? 'btn-secondary opacity-50' : 'btn-primary'
              }`}
            >
              {isCalculating ? (
                <>
                  <span className="animate-spin">⚙️</span>
                  Calculating...
                </>
              ) : (
                'Start Calculation'
              )}
            </button>
          </form>
        </div>

        {/* Results Section */}
        {results && (
          <div className="p-4 border-t border-border">
            <h3 className="font-semibold mb-3">Calculation Results</h3>
            <div className="space-y-3">
              <div className="card">
                <h4 className="font-medium mb-2">Energy Information</h4>
                <div className="space-y-1 text-sm">
                  <div className="flex justify-between">
                    <span>Method:</span>
                    <span className="font-mono">{results.method}</span>
                  </div>
                  {results.functional && (
                    <div className="flex justify-between">
                      <span>Functional:</span>
                      <span className="font-mono">{results.functional}</span>
                    </div>
                  )}
                  <div className="flex justify-between">
                    <span>Total Energy:</span>
                    <span className="font-mono">
                      {results.energy?.toFixed(8)} Ha
                    </span>
                  </div>
                  <div className="flex justify-between">
                    <span>Converged:</span>
                    <span className={results.converged ? 'text-success' : 'text-error'}>
                      {results.converged ? 'Yes' : 'No'}
                    </span>
                  </div>
                </div>
              </div>

              {results.homo_energy && (
                <div className="card">
                  <h4 className="font-medium mb-2">Orbital Information</h4>
                  <div className="space-y-1 text-sm">
                    <div className="flex justify-between">
                      <span>HOMO:</span>
                      <span className="font-mono">
                        {results.homo_energy.toFixed(4)} Ha
                      </span>
                    </div>
                    {results.lumo_energy && (
                      <>
                        <div className="flex justify-between">
                          <span>LUMO:</span>
                          <span className="font-mono">
                            {results.lumo_energy.toFixed(4)} Ha
                          </span>
                        </div>
                        <div className="flex justify-between">
                          <span>HOMO-LUMO Gap:</span>
                          <span className="font-mono">
                            {results.homo_lumo_gap.toFixed(4)} Ha
                          </span>
                        </div>
                      </>
                    )}
                  </div>
                </div>
              )}

              {results.dipole_moment && (
                <div className="card">
                  <h4 className="font-medium mb-2">Molecular Properties</h4>
                  <div className="space-y-1 text-sm">
                    <div className="flex justify-between">
                      <span>Dipole Moment:</span>
                      <span className="font-mono">
                        {results.dipole_magnitude?.toFixed(4)} Debye
                      </span>
                    </div>
                  </div>
                </div>
              )}
            </div>
          </div>
        )}
      </div>
    </div>
  )
}

export default CalculationPanel