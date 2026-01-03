import { useMutation, useQuery, useQueryClient } from '@tanstack/react-query';
import { getGpu4PyscfStatus, installGpu4Pyscf } from '../apiClient';
import { components } from '../types/generated-api';

type Gpu4PyscfStatus = components['schemas']['Gpu4PyscfStatus'];
type Gpu4PyscfInstallRequest = components['schemas']['Gpu4PyscfInstallRequest'];
type Gpu4PyscfInstallResult = components['schemas']['Gpu4PyscfInstallResult'];

const gpu4pyscfKeys = {
  all: ['gpu4pyscf'] as const,
  status: () => [...gpu4pyscfKeys.all, 'status'] as const,
};

export const useGpu4PyscfStatus = () => {
  return useQuery({
    queryKey: gpu4pyscfKeys.status(),
    queryFn: () => getGpu4PyscfStatus(),
    staleTime: 60 * 1000,
    retry: 1,
  });
};

export const useInstallGpu4Pyscf = () => {
  const queryClient = useQueryClient();

  return useMutation({
    mutationFn: (payload?: Gpu4PyscfInstallRequest) =>
      installGpu4Pyscf(payload),
    onSuccess: (result: Gpu4PyscfInstallResult) => {
      if (result?.status) {
        queryClient.setQueryData(gpu4pyscfKeys.status(), result.status);
      } else {
        queryClient.invalidateQueries({ queryKey: gpu4pyscfKeys.status() });
      }
    },
  });
};

export const useGpu4Pyscf = () => {
  const statusQuery = useGpu4PyscfStatus();
  const installMutation = useInstallGpu4Pyscf();

  return {
    status: statusQuery.data as Gpu4PyscfStatus | undefined,
    isLoading: statusQuery.isLoading,
    isFetching: statusQuery.isFetching,
    isInstalling: installMutation.isPending,
    error: statusQuery.error || installMutation.error,
    installGpu4Pyscf: installMutation.mutate,
    installGpu4PyscfAsync: installMutation.mutateAsync,
    refetch: statusQuery.refetch,
  };
};
