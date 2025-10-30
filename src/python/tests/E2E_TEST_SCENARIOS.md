# End-to-End Test Scenarios for LangGraph Multi-Agent Dispatcher

This document outlines manual E2E test scenarios to verify the LangGraph multi-agent dispatcher integration with the PySCF_front application.

## Test Environment Setup

Before running these tests:

1. Start the development server:
   ```bash
   npm run dev
   ```

2. Open the application in your browser (usually at http://localhost:3000 or the port shown in the Electron window)

3. Navigate to the Agent page (AI Assistant tab)

## Test Scenarios

### Scenario 1: Molecular Agent Routing - Quantum Chemistry Query

**Objective**: Verify that chemistry-related queries are routed to the MolecularAgent.

**Test Steps**:
1. In the Agent chat interface, enter the following query:
   ```
   水分子（H2O）のHOMO軌道エネルギーについて教えてください
   ```
   (English: "Tell me about the HOMO orbital energy of water molecule (H2O)")

2. Alternative queries to test:
   - "Calculate the energy of benzene using DFT"
   - "What is the LUMO energy of methane?"
   - "Perform a geometry optimization on H2"

**Expected Results**:
- ✅ Response should be related to quantum chemistry/molecular calculations
- ✅ Response should mention concepts like HOMO, orbitals, energy levels, or quantum chemistry
- ✅ Response should NOT attempt to search for papers (no arXiv links)
- ✅ Streaming should work smoothly (text appears progressively)
- ✅ No errors in the browser console
- ✅ Server logs should show: `"Invoking LangGraph dispatcher"` and routing to molecular agent

**Server Log Verification**:
Check the terminal for logs like:
```
[INFO] Invoking LangGraph dispatcher
[DEBUG] Received graph state update: dict_keys(['molecular'])
```

---

### Scenario 2: Research Agent Routing - Literature Search Query

**Objective**: Verify that paper search queries are routed to the ResearchAgent.

**Test Steps**:
1. In the Agent chat interface, enter the following query:
   ```
   密度汎関数理論に関する最新の論文を探してください
   ```
   (English: "Find recent papers on density functional theory")

2. Alternative queries to test:
   - "Find papers about CASSCF method"
   - "Search for literature on time-dependent DFT"
   - "What are recent publications about quantum chemistry calculations?"

**Expected Results**:
- ✅ Response should contain a list of academic papers
- ✅ Response should include paper titles, authors, publication dates
- ✅ Response should include clickable arXiv PDF links (format: `http://arxiv.org/pdf/...`)
- ✅ Response should be formatted in Markdown (headers, lists, links)
- ✅ Streaming should work smoothly
- ✅ No errors in the browser console
- ✅ Server logs should show routing to research agent

**Server Log Verification**:
Check the terminal for logs like:
```
[INFO] Invoking LangGraph dispatcher
[DEBUG] Received graph state update: dict_keys(['research'])
[DEBUG] Executing research_agent_node
```

---

### Scenario 3: Conversation History - Context Retention

**Objective**: Verify that conversation history is properly maintained and converted.

**Test Steps**:
1. Start a conversation with a molecular query:
   ```
   Tell me about water molecule
   ```

2. Follow up with a contextual query:
   ```
   What is its HOMO energy?
   ```

3. Then switch to research:
   ```
   Find papers about it
   ```

**Expected Results**:
- ✅ Second message should understand "its" refers to water
- ✅ Third message should understand "it" refers to water molecule
- ✅ Responses should be contextually appropriate
- ✅ No loss of conversation context

**Server Log Verification**:
```
[DEBUG] Converted X history messages to LangChain format
```

---

### Scenario 4: Error Handling - Invalid Query

**Objective**: Verify graceful error handling.

**Test Steps**:
1. Enter an empty message (should be prevented by UI)
2. Enter an extremely long message (>100,000 characters)
3. Disconnect network mid-stream (if testing network resilience)

**Expected Results**:
- ✅ Empty messages should be blocked by frontend validation
- ✅ Long messages should return a 400 error with appropriate message
- ✅ Network errors should display user-friendly error message
- ✅ No application crashes

---

### Scenario 5: Streaming Behavior

**Objective**: Verify SSE streaming works correctly.

**Test Steps**:
1. Enter any query that generates a long response
2. Observe the response as it arrives

**Expected Results**:
- ✅ Text should appear progressively (not all at once)
- ✅ "Done" event should be received at the end
- ✅ No duplicate text or garbled characters
- ✅ Markdown formatting should render correctly

**Browser Console Verification**:
Open Developer Tools > Console, check for:
```
[SSE Debug] Starting SSE stream
[SSE Debug] Received SSE message { type: 'chunk', hasText: true }
[SSE Debug] Stream completed
```

---

### Scenario 6: Router Decision Validation

**Objective**: Verify router makes correct decisions for ambiguous queries.

**Test Steps**:
Test borderline queries that could go either way:

1. "What is DFT?" (Should go to Molecular - explanation)
2. "Who invented DFT?" (Could go to Research - historical)
3. "Calculate DFT for H2" (Should go to Molecular - calculation)
4. "Papers on DFT applications" (Should go to Research - literature)

**Expected Results**:
- ✅ Router should make reasonable decisions
- ✅ Responses should be appropriate for the chosen agent
- ✅ No routing errors

---

## Integration Test Verification

After manual testing, run the automated integration tests:

```bash
cd src/python
conda activate pyscf-env
pytest tests/integration/test_api_endpoints/test_agent_api.py -v
```

**Expected Output**:
```
tests/integration/test_api_endpoints/test_agent_api.py::TestAgentChatAPI::test_chat_success_with_streaming PASSED
tests/integration/test_api_endpoints/test_agent_api.py::TestAgentChatAPI::test_chat_with_history PASSED
tests/integration/test_api_endpoints/test_agent_api.py::TestAgentChatAPI::test_chat_empty_message PASSED
tests/integration/test_api_endpoints/test_agent_api.py::TestAgentChatAPI::test_chat_routing_to_molecular_agent PASSED
tests/integration/test_api_endpoints/test_agent_api.py::TestAgentChatAPI::test_chat_routing_to_research_agent PASSED
...
```

All tests should PASS.

---

## Success Criteria

The Phase 4 implementation is considered successful when:

- [ ] All integration tests pass
- [ ] Molecular queries route to MolecularAgent and return chemistry-related responses
- [ ] Research queries route to ResearchAgent and return formatted paper listings with arXiv links
- [ ] Conversation history is properly maintained across messages
- [ ] SSE streaming works without errors
- [ ] Error handling is graceful and user-friendly
- [ ] No regression in existing functionality
- [ ] Server logs show proper graph execution flow

---

## Troubleshooting

### Common Issues

**Issue**: "Graph execution error" in response
- **Solution**: Check that all dependencies (langchain, langgraph, arxiv) are installed
- **Verify**: Run `conda list | grep langchain`

**Issue**: Routing always goes to the same agent
- **Solution**: Check router_node implementation in `src/python/agent/graph.py`
- **Verify**: Add debug logging to see router decision

**Issue**: Streaming stops mid-response
- **Solution**: Check for exceptions in server logs
- **Verify**: Look for Python tracebacks in terminal output

**Issue**: "No AI response extracted from graph"
- **Solution**: Check that agent nodes are properly appending AIMessage to state
- **Verify**: Add logging in `_create_graph_stream` to inspect final_state

---

## Next Steps

After successful E2E testing:

1. Document any edge cases discovered during testing
2. Add unit tests for any new edge cases
3. Update user documentation if new features are exposed
4. Consider adding telemetry/analytics for router decisions (for future improvements)
