type MainRendererEntryParams = {
  backendPort: number;
  htmlPath: string;
  isPackaged: boolean;
  rendererUrl?: string;
};

type SplashRendererEntryParams = {
  htmlPath: string;
  isPackaged: boolean;
  rendererUrl?: string;
};

type FileRendererEntry = {
  type: 'file';
  path: string;
  query?: Record<string, string>;
};

type UrlRendererEntry = {
  type: 'url';
  url: string;
};

export type RendererEntry = FileRendererEntry | UrlRendererEntry;

const buildDevUrl = (rendererUrl: string, pathname: string): URL => {
  const baseUrl = rendererUrl.endsWith('/') ? rendererUrl : `${rendererUrl}/`;
  return new URL(pathname, baseUrl);
};

export const getMainRendererEntry = ({
  backendPort,
  htmlPath,
  isPackaged,
  rendererUrl,
}: MainRendererEntryParams): RendererEntry => {
  if (!isPackaged && rendererUrl) {
    const url = buildDevUrl(rendererUrl, '/');
    url.searchParams.set('backend_port', String(backendPort));
    return {
      type: 'url',
      url: url.toString(),
    };
  }

  return {
    type: 'file',
    path: htmlPath,
    query: {
      backend_port: String(backendPort),
    },
  };
};

export const getSplashRendererEntry = ({
  htmlPath,
  isPackaged,
  rendererUrl,
}: SplashRendererEntryParams): RendererEntry => {
  if (!isPackaged && rendererUrl) {
    return {
      type: 'url',
      url: buildDevUrl(rendererUrl, 'splash.html').toString(),
    };
  }

  return {
    type: 'file',
    path: htmlPath,
  };
};
