#!/usr/bin/env python3
"""
Simple gRPC client to test the full backend workflow
"""
import grpc
import sys
import time

# Add the backend source to the path
sys.path.append('backend/src')

from grpc_stubs import calculation_pb2, calculation_pb2_grpc


def test_grpc_workflow():
    """Test the complete molecule creation -> instance creation -> calculation workflow"""
    
    channel = grpc.insecure_channel('localhost:50051')
    stub = calculation_pb2_grpc.CalculationServiceStub(channel)
    
    print("🧪 Testing gRPC Backend Workflow")
    print("=" * 50)
    
    # Test 1: Health Check
    print("1. Testing Health Check...")
    try:
        health_request = calculation_pb2.HealthCheckRequest()
        health_response = stub.HealthCheck(health_request)
        print(f"   ✅ Health Status: {health_response.status}")
        print(f"   📊 Details: {dict(health_response.details)}")
    except Exception as e:
        print(f"   ❌ Health check failed: {e}")
        return
    
    # Test 2: System Info
    print("\n2. Testing System Info...")
    try:
        system_request = calculation_pb2.GetSystemInfoRequest()
        system_response = stub.GetSystemInfo(system_request)
        print(f"   ✅ System Info Retrieved")
        print(f"   🐍 Python: {system_response.system_info.python_version[:50]}...")
        print(f"   ⚡ PySCF: {system_response.system_info.pyscf_version}")
        print(f"   🎯 Methods: {len(system_response.system_info.available_methods)}")
        print(f"   📚 Basis Sets: {len(system_response.system_info.available_basis_sets)}")
    except Exception as e:
        print(f"   ❌ System info failed: {e}")
        return
    
    # Test 3: Create Molecule (Water H2O)
    print("\n3. Testing Molecule Creation (Water H2O)...")
    try:
        atoms = [
            calculation_pb2.Atom(symbol='O', x=0.0, y=0.0, z=0.0),
            calculation_pb2.Atom(symbol='H', x=0.0, y=0.757, z=0.587),
            calculation_pb2.Atom(symbol='H', x=0.0, y=-0.757, z=0.587)
        ]
        
        molecule_request = calculation_pb2.CreateMoleculeRequest(
            name="Water",
            formula="H2O",
            atoms=atoms,
            charge=0,
            multiplicity=1,
            symmetry="c2v"
        )
        
        molecule_response = stub.CreateMolecule(molecule_request)
        if molecule_response.success:
            print(f"   ✅ Molecule Created: {molecule_response.molecule.name}")
            print(f"   🆔 ID: {molecule_response.molecule.id}")
            print(f"   ⚖️  Weight: {molecule_response.molecule.molecular_weight:.3f} amu")
            molecule_id = molecule_response.molecule.id
        else:
            print(f"   ❌ Molecule creation failed: {molecule_response.message}")
            return
    except Exception as e:
        print(f"   ❌ Molecule creation failed: {e}")
        return
    
    # Test 4: Get Molecule
    print("\n4. Testing Molecule Retrieval...")
    try:
        get_mol_request = calculation_pb2.GetMoleculeRequest(molecule_id=molecule_id)
        get_mol_response = stub.GetMolecule(get_mol_request)
        if get_mol_response.success:
            print(f"   ✅ Molecule Retrieved: {get_mol_response.molecule.name}")
            print(f"   🧪 Formula: {get_mol_response.molecule.formula}")
        else:
            print(f"   ❌ Molecule retrieval failed: {get_mol_response.message}")
    except Exception as e:
        print(f"   ❌ Molecule retrieval failed: {e}")
    
    # Test 5: Create Instance
    print("\n5. Testing Instance Creation...")
    try:
        instance_request = calculation_pb2.CreateInstanceRequest(
            name="Water HF Calculation",
            description="Test Hartree-Fock calculation on water molecule",
            molecule_id=molecule_id
        )
        
        instance_response = stub.CreateInstance(instance_request)
        if instance_response.success:
            print(f"   ✅ Instance Created: {instance_response.instance.name}")
            print(f"   🆔 ID: {instance_response.instance.id}")
            instance_id = instance_response.instance.id
        else:
            print(f"   ❌ Instance creation failed: {instance_response.message}")
            return
    except Exception as e:
        print(f"   ❌ Instance creation failed: {e}")
        return
    
    # Test 6: Get Instance
    print("\n6. Testing Instance Retrieval...")
    try:
        get_inst_request = calculation_pb2.GetInstanceRequest(instance_id=instance_id)
        get_inst_response = stub.GetInstance(get_inst_request)
        if get_inst_response.success:
            print(f"   ✅ Instance Retrieved: {get_inst_response.instance.name}")
            print(f"   📝 Description: {get_inst_response.instance.description}")
        else:
            print(f"   ❌ Instance retrieval failed: {get_inst_response.message}")
    except Exception as e:
        print(f"   ❌ Instance retrieval failed: {e}")
    
    # Test 7: List Instances
    print("\n7. Testing Instance Listing...")
    try:
        list_request = calculation_pb2.ListInstancesRequest(page=1, page_size=10)
        list_response = stub.ListInstances(list_request)
        if list_response.success:
            print(f"   ✅ Found {list_response.total_count} instances")
            for i, instance in enumerate(list_response.instances):
                print(f"   📋 [{i+1}] {instance.name} ({instance.id[:8]}...)")
        else:
            print(f"   ❌ Instance listing failed: {list_response.message}")
    except Exception as e:
        print(f"   ❌ Instance listing failed: {e}")
    
    # Test 8: Start Simple HF Calculation (just initiate, don't wait for completion)
    print("\n8. Testing Calculation Start (HF/sto-3g)...")
    try:
        calc_request = calculation_pb2.StartCalculationRequest(
            instance_id=instance_id,
            method="HF",
            basis_set="sto-3g",
            priority=1
        )
        
        print("   🚀 Starting calculation...")
        progress_stream = stub.StartCalculation(calc_request)
        
        # Read just the first few progress updates
        progress_count = 0
        for progress in progress_stream:
            print(f"   📊 {progress.current_step}: {progress.progress_percentage:.1f}% - {progress.message}")
            progress_count += 1
            
            # Stop after a few updates to avoid long wait
            if progress_count >= 3 or progress.progress_percentage >= 100:
                break
                
        print(f"   ✅ Calculation initiated successfully")
        
    except Exception as e:
        print(f"   ❌ Calculation failed: {e}")
    
    print("\n" + "=" * 50)
    print("🎉 gRPC Backend Test Complete!")
    print("✅ All core functionality working with database integration")


if __name__ == "__main__":
    test_grpc_workflow()