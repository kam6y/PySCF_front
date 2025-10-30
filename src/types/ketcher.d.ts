// Type definitions for ketcher packages
declare module 'ketcher-standalone' {
  export class StandaloneStructServiceProvider {
    constructor();
  }
}

declare module 'ketcher-react' {
  import { ReactElement } from 'react';
  import { Ketcher } from 'ketcher-core';

  export interface EditorProps {
    staticResourcesUrl?: string;
    structServiceProvider: any;
    errorHandler?: (message: string) => void;
    onInit?: (ketcher: Ketcher) => void;
  }

  export const Editor: (props: EditorProps) => ReactElement;
}

declare module 'ketcher-core' {
  export class Ketcher {
    getSmiles(isExtended?: boolean): Promise<string>;
    getMolfile(molfileFormat?: 'v2000' | 'v3000' | 'auto'): Promise<string>;
    getRxn(molfileFormat?: 'v2000' | 'v3000'): Promise<string>;
    setMolecule(molecule: string | object): Promise<void>;
    editor: {
      subscribe: (event: string, callback: (operations: any) => void) => void;
    };
  }
}

declare module 'miew' {
  const Miew: any;
  export default Miew;
}
