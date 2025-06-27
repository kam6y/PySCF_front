import React, { useState } from 'react'

interface ProjectPanelProps {
  onMoleculeLoad: (moleculeData: any) => void
  currentMolecule: any
}

const ProjectPanel: React.FC<ProjectPanelProps> = ({ 
  onMoleculeLoad, 
  currentMolecule 
}) => {
  const [activeTab, setActiveTab] = useState('molecules')

  // Sample molecules for testing
  const sampleMolecules = [
    {
      name: 'Water (H₂O)',
      formula: 'H₂O',
      data: {
        type: 'coordinates',
        coordinates: [
          ['O', 0.0, 0.0, 0.0],
          ['H', 0.757, 0.586, 0.0],
          ['H', -0.757, 0.586, 0.0]
        ],
        charge: 0,
        spin: 0
      }
    },
    {
      name: 'Methane (CH₄)',
      formula: 'CH₄',
      data: {
        type: 'coordinates',
        coordinates: [
          ['C', 0.0, 0.0, 0.0],
          ['H', 1.089, 0.0, 0.0],
          ['H', -0.363, 1.027, 0.0],
          ['H', -0.363, -0.513, 0.889],
          ['H', -0.363, -0.513, -0.889]
        ],
        charge: 0,
        spin: 0
      }
    },
    {
      name: 'Hydrogen (H₂)',
      formula: 'H₂',
      data: {
        type: 'coordinates',
        coordinates: [
          ['H', 0.0, 0.0, 0.0],
          ['H', 0.74, 0.0, 0.0]
        ],
        charge: 0,
        spin: 0
      }
    }
  ]

  const handleMoleculeSelect = (molecule: any) => {
    onMoleculeLoad(molecule.data)
  }

  return (
    <div className="h-full flex flex-col">
      {/* Header */}
      <div className="p-4 border-b border-border">
        <h2 className="text-lg font-semibold">PySCF Front</h2>
        <p className="text-sm text-secondary">Quantum Chemistry GUI</p>
      </div>

      {/* Tab Navigation */}
      <div className="flex border-b border-border">
        <button
          onClick={() => setActiveTab('molecules')}
          className={`flex-1 py-3 px-4 text-sm font-medium ${
            activeTab === 'molecules'
              ? 'border-b-2 border-primary text-primary'
              : 'text-secondary hover:text-primary'
          }`}
        >
          Molecules
        </button>
        <button
          onClick={() => setActiveTab('projects')}
          className={`flex-1 py-3 px-4 text-sm font-medium ${
            activeTab === 'projects'
              ? 'border-b-2 border-primary text-primary'
              : 'text-secondary hover:text-primary'
          }`}
        >
          Projects
        </button>
      </div>

      <div className="flex-1 overflow-y-auto">
        {activeTab === 'molecules' && (
          <div className="p-4">
            <div className="mb-4">
              <h3 className="font-medium mb-2">Sample Molecules</h3>
              <p className="text-sm text-secondary mb-4">
                Select a molecule to load into the viewer
              </p>
            </div>

            <div className="space-y-2">
              {sampleMolecules.map((molecule, index) => (
                <div
                  key={index}
                  onClick={() => handleMoleculeSelect(molecule)}
                  className={`p-3 rounded-lg border cursor-pointer transition-colors ${
                    currentMolecule && 
                    JSON.stringify(currentMolecule.coordinates) === 
                    JSON.stringify(molecule.data.coordinates)
                      ? 'border-primary bg-primary bg-opacity-10'
                      : 'border-border hover:border-primary hover:bg-tertiary'
                  }`}
                >
                  <div className="flex items-center justify-between">
                    <div>
                      <h4 className="font-medium">{molecule.name}</h4>
                      <p className="text-sm text-secondary">{molecule.formula}</p>
                    </div>
                    <div className="text-2xl">⚛️</div>
                  </div>
                </div>
              ))}
            </div>

            <div className="mt-6">
              <h3 className="font-medium mb-2">Input Methods</h3>
              <div className="space-y-2">
                <button className="btn btn-secondary w-full text-sm">
                  📁 Load XYZ File
                </button>
                <button className="btn btn-secondary w-full text-sm">
                  🧪 Draw Molecule
                </button>
                <button className="btn btn-secondary w-full text-sm">
                  📝 SMILES Input
                </button>
                <button className="btn btn-secondary w-full text-sm">
                  🔍 PubChem Search
                </button>
              </div>
            </div>
          </div>
        )}

        {activeTab === 'projects' && (
          <div className="p-4">
            <div className="mb-4">
              <h3 className="font-medium mb-2">Recent Projects</h3>
              <p className="text-sm text-secondary mb-4">
                Manage your calculation projects
              </p>
            </div>

            <div className="text-center text-muted py-8">
              <div className="text-4xl mb-2">📁</div>
              <p>No projects yet</p>
              <p className="text-sm">
                Run calculations to create projects
              </p>
            </div>

            <div className="space-y-2">
              <button className="btn btn-primary w-full text-sm">
                ➕ New Project
              </button>
              <button className="btn btn-secondary w-full text-sm">
                📂 Open Project
              </button>
            </div>
          </div>
        )}
      </div>
    </div>
  )
}

export default ProjectPanel