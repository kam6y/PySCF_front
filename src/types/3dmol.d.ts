export interface AtomSpec {
  elem?: string;
  chain?: string;
  resi?: number;
  resn?: string;
  atom?: string;
  serial?: number;
  altLoc?: string;
  x?: number;
  y?: number;
  z?: number;
  properties?: {
    [key: string]: any;
    charge?: number;
  };
}

export interface StyleSpec {
  line?: {
    linewidth?: number;
    color?: string;
  };
  stick?: {
    radius?: number;
    color?: string;
    colorscheme?: string;
  };
  sphere?: {
    radius?: number;
    color?: string;
    colorscheme?: string;
  };
  cartoon?: {
    color?: string;
    colorscheme?: string;
    style?: string;
  };
  surface?: {
    opacity?: number;
    color?: string;
    colorscheme?: string;
  };
}

export interface ViewerSpec {
  defaultcolors?: any;
  nomouse?: boolean;
  backgroundColor?: string;
}

export interface LabelSpec {
  position?: { x: number; y: number; z: number } | { serial: number };
  fontColor?: string;
  backgroundColor?: string;
  backgroundOpacity?: number;
  showBackground?: boolean;
  fontSize?: number;
  inFront?: boolean;
  bold?: boolean;
}

export interface Label {
  setText(text: string): void;
  remove(): void;
}

export interface GLViewer {
  addModel(data: string, format?: string, options?: any): GLModel;
  removeModel(model: GLModel): void;
  removeAllModels(): void;
  createModelFrom(sel: AtomSpec, extract?: boolean): GLModel;

  setStyle(sel: AtomSpec, style: StyleSpec): GLViewer;
  addStyle(sel: AtomSpec, style: StyleSpec): GLViewer;

  zoomTo(
    sel?: AtomSpec,
    animationDuration?: number,
    fixedPath?: boolean
  ): GLViewer;
  zoom(factor?: number, animationDuration?: number): GLViewer;
  translate(x: number, y: number, animationDuration?: number): GLViewer;
  rotate(angle: number, axis?: string, animationDuration?: number): GLViewer;

  render(callback?: () => void): GLViewer;
  clear(): GLViewer;

  setBackgroundColor(color: string, alpha?: number): GLViewer;
  setWidth(width: number): GLViewer;
  setHeight(height: number): GLViewer;
  resize(): GLViewer;

  getView(): number[];
  setView(view: number[]): GLViewer;

  addSurface(
    type: string,
    style: any,
    sel?: AtomSpec & { charges?: number[] },
    allsel?: AtomSpec
  ): number;
  removeSurface(surfid: number): GLViewer;
  removeAllSurfaces(): GLViewer;

  spin(axis?: string, speed?: number): GLViewer;
  stopAnimate(): GLViewer;

  screenshot(width?: number, height?: number, format?: string): string;
  getModel(id?: number): GLModel;

  addArrow(spec: any): void;
  removeAllShapes(): GLViewer;
  addLabel(text: string, options: LabelSpec): Label;
  removeAllLabels(): GLViewer;
}

export interface GLModel {
  setStyle(sel: AtomSpec, style: StyleSpec, add?: boolean): GLModel;
  addStyle(sel: AtomSpec, style: StyleSpec): GLModel;
  removeStyle(sel: AtomSpec, style: StyleSpec): GLModel;

  selectedAtoms(sel: AtomSpec): AtomSpec[];
  atoms: AtomSpec[];

  setColorByElement(sel: AtomSpec, colors?: any): GLModel;
  setColorByProperty(sel: AtomSpec, prop: string, scheme?: any): GLModel;

  addSurface(
    type: string,
    style: any,
    sel?: AtomSpec & { charges?: number[] },
    allsel?: AtomSpec
  ): number;
  removeSurface(surf: number): GLModel;

  addLine(spec: any): GLModel;
  addArrow(spec: any): GLModel;
  addSphere(spec: any): GLModel;
  addCylinder(spec: any): GLModel;
  addCurve(spec: any): GLModel;
}

declare module '3dmol' {
  export function createViewer(
    element: HTMLElement | string,
    config?: ViewerSpec
  ): GLViewer;
  export function createViewerGrid(
    element: HTMLElement | string,
    config?: ViewerSpec,
    rows?: number,
    cols?: number
  ): GLViewer[];

  export function download(
    query: string,
    viewer: GLViewer,
    options?: any,
    callback?: () => void
  ): void;

  export const ElementColors: { [element: string]: number };

  export enum SurfaceType {
    VDW = 1,
    MS = 2,
    SAS = 3,
    SES = 4,
  }

  export namespace Gradient {
    class RWB {
      constructor(min: number, max: number);
    }
    class ROYGB {
      constructor(min: number, max: number);
    }
    class Sinebow {
      constructor(min: number, max: number);
    }
  }
}
