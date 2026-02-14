import {
  PubChemSearchResponseData,
  SMILESConvertResponseData,
} from '../types/api-types';
import { request } from './core';

export const searchPubChem = (
  query: string,
  searchType: 'name' | 'cid'
): Promise<PubChemSearchResponseData> => {
  return request<PubChemSearchResponseData>('/api/pubchem/search', {
    method: 'POST',
    body: JSON.stringify({ query, searchType }),
  });
};

export const convertSmilesToXyz = (
  smiles: string
): Promise<SMILESConvertResponseData> => {
  return request<SMILESConvertResponseData>('/api/smiles/convert', {
    method: 'POST',
    body: JSON.stringify({ smiles }),
  });
};
