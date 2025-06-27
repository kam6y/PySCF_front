"""
Integration tests for Python backend stdin/stdout communication
Tests the main communication loop and JSON message processing
"""

import unittest  
import sys
import os
import json
import subprocess
import threading
import time
import io
from unittest.mock import patch, MagicMock
from io import StringIO

# Add src/python to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/python'))

from main import PySCFBackend


class TestPythonStdinCommunication(unittest.TestCase):
    """Test Python backend stdin/stdout communication protocol"""
    
    def setUp(self):
        """Set up test fixtures"""
        self.backend = PySCFBackend()
        
        # Test messages
        self.test_messages = {
            'valid_status': {
                'action': 'get-status',
                'data': {},
                'id': 'test_1'
            },
            'valid_calculate': {
                'action': 'calculate',
                'data': {
                    'molecule': {
                        'type': 'coordinates',
                        'coordinates': [
                            ['H', 0.0, 0.0, 0.0],
                            ['H', 0.74, 0.0, 0.0]
                        ],
                        'charge': 0,
                        'spin': 0
                    },
                    'method': 'HF',
                    'basis': 'STO-3G'
                },
                'id': 'test_calc'
            },
            'invalid_json': '{"invalid": json}',
            'missing_action': {
                'data': {},
                'id': 'test_2'
            },
            'unknown_action': {
                'action': 'unknown_command',
                'data': {},
                'id': 'test_3'
            }
        }
    
    def test_handle_valid_status_message(self):
        """Test handling valid status message"""
        message = self.test_messages['valid_status']
        response = self.backend.handle_message(message)
        
        self.assertEqual(response['id'], 'test_1')
        self.assertEqual(response['action'], 'get-status')
        self.assertEqual(response['status'], 'success')
        self.assertEqual(response['backend_status'], 'running')
        self.assertIn('jobs', response)
    
    def test_handle_valid_calculation_message(self):
        """Test handling valid calculation message"""
        message = self.test_messages['valid_calculate']
        response = self.backend.handle_message(message)
        
        self.assertEqual(response['id'], 'test_calc')
        self.assertEqual(response['action'], 'calculate')
        self.assertEqual(response['status'], 'success')
        self.assertIn('results', response)
        
        results = response['results']
        self.assertEqual(results['method'], 'HF')
        self.assertTrue(results['converged'])
        self.assertLess(results['energy'], 0.0)
    
    def test_handle_missing_action_message(self):
        """Test handling message with missing action field"""
        message = self.test_messages['missing_action']
        response = self.backend.handle_message(message)
        
        self.assertEqual(response['id'], 'test_2')
        self.assertEqual(response['status'], 'error')
        self.assertIn('message', response)
    
    def test_handle_unknown_action_message(self):
        """Test handling message with unknown action"""
        message = self.test_messages['unknown_action']
        response = self.backend.handle_message(message)
        
        self.assertEqual(response['id'], 'test_3')
        self.assertEqual(response['status'], 'error')
        self.assertIn('Unknown action', response['message'])
    
    def test_handle_malformed_message(self):
        """Test handling completely malformed message"""
        malformed_message = {'invalid_structure': True}
        response = self.backend.handle_message(malformed_message)
        
        self.assertEqual(response['status'], 'error')
        self.assertIn('message', response)
    
    def test_stop_message_sets_running_false(self):
        """Test that stop message sets running flag to False"""
        message = {
            'action': 'stop',
            'data': {},
            'id': 'stop_test'
        }
        
        self.assertTrue(self.backend.running)
        response = self.backend.handle_message(message)
        
        self.assertFalse(self.backend.running)
        self.assertEqual(response['status'], 'success')
        self.assertEqual(response['message'], 'Backend stopping')


class TestPythonMainLoop(unittest.TestCase):
    """Test the main stdin processing loop"""
    
    def setUp(self):
        """Set up test fixtures"""
        self.backend = PySCFBackend()
    
    @patch('sys.stdin')
    @patch('builtins.print')
    def test_main_loop_processes_single_message(self, mock_print, mock_stdin):
        """Test main loop processes single message correctly"""
        # Mock stdin to return a single message then EOF
        test_message = {
            'action': 'get-status',
            'data': {},
            'id': 'loop_test_1'
        }
        
        mock_stdin.readline.side_effect = [
            json.dumps(test_message) + '\n',
            ''  # EOF
        ]
        
        # Run the main loop
        self.backend.run()
        
        # Verify print was called with JSON response
        self.assertEqual(mock_print.call_count, 1)
        printed_output = mock_print.call_args[0][0]
        response = json.loads(printed_output)
        
        self.assertEqual(response['id'], 'loop_test_1')
        self.assertEqual(response['status'], 'success')
    
    @patch('sys.stdin')
    @patch('builtins.print')
    def test_main_loop_processes_multiple_messages(self, mock_print, mock_stdin):
        """Test main loop processes multiple messages in sequence"""
        # Mock stdin to return multiple messages
        messages = [
            {'action': 'get-status', 'data': {}, 'id': 'msg_1'},
            {'action': 'get-status', 'data': {}, 'id': 'msg_2'},
            {'action': 'stop', 'data': {}, 'id': 'msg_3'}
        ]
        
        mock_stdin.readline.side_effect = [
            json.dumps(msg) + '\n' for msg in messages
        ] + ['']  # EOF
        
        # Run the main loop
        self.backend.run()
        
        # Verify all messages were processed
        self.assertEqual(mock_print.call_count, 3)
        
        # Check responses
        for i, call in enumerate(mock_print.call_args_list):
            printed_output = call[0][0]
            response = json.loads(printed_output)
            self.assertEqual(response['id'], f'msg_{i+1}')
    
    @patch('sys.stdin')
    @patch('builtins.print')
    def test_main_loop_handles_json_decode_error(self, mock_print, mock_stdin):
        """Test main loop handles JSON decode errors gracefully"""
        # Mock stdin to return invalid JSON
        mock_stdin.readline.side_effect = [
            'invalid json string\n',
            ''  # EOF
        ]
        
        # Run the main loop
        self.backend.run()
        
        # Verify error response was printed
        self.assertEqual(mock_print.call_count, 1)
        printed_output = mock_print.call_args[0][0]
        response = json.loads(printed_output)
        
        self.assertEqual(response['status'], 'error')
        self.assertEqual(response['message'], 'Invalid JSON format')
    
    @patch('sys.stdin')
    @patch('builtins.print')
    def test_main_loop_stops_on_stop_message(self, mock_print, mock_stdin):
        """Test main loop stops when receiving stop message"""
        # Mock stdin to return stop message
        stop_message = {
            'action': 'stop',
            'data': {},
            'id': 'stop_test'
        }
        
        mock_stdin.readline.side_effect = [
            json.dumps(stop_message) + '\n',
            json.dumps({'action': 'get-status', 'data': {}, 'id': 'after_stop'}) + '\n',
            ''  # EOF
        ]
        
        # Run the main loop
        self.backend.run()
        
        # Should only process the stop message, not the subsequent one
        self.assertEqual(mock_print.call_count, 1)
        printed_output = mock_print.call_args[0][0]
        response = json.loads(printed_output)
        
        self.assertEqual(response['id'], 'stop_test')
        self.assertEqual(response['status'], 'success')
    
    @patch('sys.stdin')
    @patch('builtins.print')
    @patch('logging.Logger.error')
    def test_main_loop_handles_unexpected_errors(self, mock_logger, mock_print, mock_stdin):
        """Test main loop handles unexpected errors in message processing"""
        # Mock stdin to return message that will cause error
        test_message = {
            'action': 'get-status',
            'data': {},
            'id': 'error_test'
        }
        
        mock_stdin.readline.side_effect = [
            json.dumps(test_message) + '\n',
            ''  # EOF
        ]
        
        # Mock handle_message to raise exception
        original_handle = self.backend.handle_message
        def error_handle(msg):
            if msg.get('id') == 'error_test':
                raise RuntimeError("Test error")
            return original_handle(msg)
        
        self.backend.handle_message = error_handle
        
        # Run the main loop
        self.backend.run()
        
        # Verify error response was printed
        self.assertEqual(mock_print.call_count, 1)
        printed_output = mock_print.call_args[0][0]
        response = json.loads(printed_output)
        
        self.assertEqual(response['status'], 'error')
        # After simplification, the error message is now generic
        self.assertEqual(response['message'], 'Message processing error')


class TestPythonSubprocessCommunication(unittest.TestCase):
    """Test communication with Python backend as subprocess"""
    
    def setUp(self):
        """Set up subprocess test"""
        self.python_executable = sys.executable
        self.backend_script = os.path.join(
            os.path.dirname(__file__),
            '../../src/python/main.py'
        )
    
    def test_subprocess_startup_and_shutdown(self):
        """Test that Python subprocess starts and shuts down properly"""
        process = subprocess.Popen(
            [self.python_executable, self.backend_script],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        
        try:
            # Send stop message
            stop_message = {
                'action': 'stop',
                'data': {},
                'id': 'shutdown_test'
            }
            
            message_str = json.dumps(stop_message) + '\n'
            process.stdin.write(message_str)
            process.stdin.flush()
            
            # Read response
            response_line = process.stdout.readline()
            response = json.loads(response_line)
            
            self.assertEqual(response['id'], 'shutdown_test')
            self.assertEqual(response['status'], 'success')
            
            # Process should exit after stop message
            return_code = process.wait(timeout=5)
            self.assertEqual(return_code, 0)
            
        except subprocess.TimeoutExpired:
            process.kill()
            self.fail("Process did not exit after stop message")
        except Exception as e:
            process.kill()  
            raise e
    
    def test_subprocess_multiple_message_exchange(self):
        """Test multiple message exchanges with subprocess"""
        process = subprocess.Popen(
            [self.python_executable, self.backend_script],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        
        try:
            messages = [
                {'action': 'get-status', 'data': {}, 'id': 'msg_1'},
                {'action': 'get-status', 'data': {}, 'id': 'msg_2'},
                {'action': 'stop', 'data': {}, 'id': 'msg_3'}
            ]
            
            responses = []
            
            for message in messages:
                message_str = json.dumps(message) + '\n'
                process.stdin.write(message_str)
                process.stdin.flush()
                
                response_line = process.stdout.readline()
                response = json.loads(response_line)
                responses.append(response)
            
            # Verify all responses
            for i, response in enumerate(responses):
                self.assertEqual(response['id'], f'msg_{i+1}')
                self.assertEqual(response['status'], 'success')
            
            # Wait for process to exit
            return_code = process.wait(timeout=5)
            self.assertEqual(return_code, 0)
            
        except subprocess.TimeoutExpired:
            process.kill()
            self.fail("Process communication failed")
        except Exception as e:
            process.kill()
            raise e
    
    def test_subprocess_calculation_workflow(self):
        """Test full calculation workflow through subprocess"""
        process = subprocess.Popen(
            [self.python_executable, self.backend_script],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        
        try:
            # Send calculation message
            calc_message = {
                'action': 'calculate',
                'data': {
                    'molecule': {
                        'type': 'coordinates',
                        'coordinates': [
                            ['H', 0.0, 0.0, 0.0],
                            ['H', 0.74, 0.0, 0.0]
                        ],
                        'charge': 0,
                        'spin': 0
                    },
                    'method': 'HF',
                    'basis': 'STO-3G'
                },
                'id': 'subprocess_calc'
            }
            
            message_str = json.dumps(calc_message) + '\n'
            process.stdin.write(message_str)
            process.stdin.flush()
            
            # Read calculation response (may take time)
            response_line = process.stdout.readline()
            response = json.loads(response_line)
            
            self.assertEqual(response['id'], 'subprocess_calc')
            self.assertEqual(response['status'], 'success')
            self.assertIn('results', response)
            
            results = response['results']
            self.assertEqual(results['method'], 'HF')
            self.assertTrue(results['converged'])
            self.assertLess(results['energy'], 0.0)
            
            # Send stop message
            stop_message = {
                'action': 'stop',
                'data': {},
                'id': 'stop_after_calc'
            }
            
            message_str = json.dumps(stop_message) + '\n'
            process.stdin.write(message_str)
            process.stdin.flush()
            
            # Read stop response
            response_line = process.stdout.readline()
            response = json.loads(response_line)
            
            self.assertEqual(response['status'], 'success')
            
            # Wait for process to exit
            return_code = process.wait(timeout=10)
            self.assertEqual(return_code, 0)
            
        except subprocess.TimeoutExpired:
            process.kill()
            self.fail("Calculation subprocess communication failed")
        except Exception as e:
            process.kill()
            raise e
    
    def test_subprocess_handles_malformed_input(self):
        """Test subprocess handles malformed input gracefully"""
        process = subprocess.Popen(
            [self.python_executable, self.backend_script],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        
        try:
            # Send malformed JSON
            malformed_input = "invalid json string\n"
            process.stdin.write(malformed_input)
            process.stdin.flush()
            
            # Should get error response
            response_line = process.stdout.readline()
            response = json.loads(response_line)
            
            self.assertEqual(response['status'], 'error')
            self.assertEqual(response['message'], 'Invalid JSON format')
            
            # Send stop message to clean shutdown
            stop_message = json.dumps({'action': 'stop', 'data': {}, 'id': 'cleanup'}) + '\n'
            process.stdin.write(stop_message)
            process.stdin.flush()
            
            # Wait for process to exit
            return_code = process.wait(timeout=5)
            self.assertEqual(return_code, 0)
            
        except subprocess.TimeoutExpired:
            process.kill()
            self.fail("Subprocess did not handle malformed input properly")
        except Exception as e:
            process.kill()
            raise e


if __name__ == '__main__':
    unittest.main()