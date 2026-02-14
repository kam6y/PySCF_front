import React from 'react';
import { ConfirmationModal } from '../ConfirmationModal';
import type { DeleteModalsState } from './useDeleteModals';

export type DeleteModalsProps = DeleteModalsState;

export const DeleteModals: React.FC<DeleteModalsProps> = ({
  isDeleteModalOpen,
  calculationToDelete,
  handleConfirmDelete,
  handleCancelDelete,
  isBulkDeleteModalOpen,
  calculationsToDelete,
  handleConfirmBulkDelete,
  handleCancelBulkDelete,
  isChatDeleteModalOpen,
  chatToDelete,
  handleConfirmChatDelete,
  handleCancelChatDelete,
  deleteChatSessionPending,
}) => {
  return (
    <>
      <ConfirmationModal
        isOpen={isDeleteModalOpen}
        title="Delete Calculation"
        message={`Are you sure you want to delete "${calculationToDelete?.name}"? This action cannot be undone.`}
        confirmButtonText="Delete"
        cancelButtonText="Cancel"
        onConfirm={handleConfirmDelete}
        onCancel={handleCancelDelete}
        isLoading={false}
      />

      <ConfirmationModal
        isOpen={isBulkDeleteModalOpen}
        title="Bulk Delete Error Calculations"
        message={`Do you want to delete ${calculationsToDelete.length} error calculations? This operation cannot be undone.`}
        confirmButtonText="Delete All"
        cancelButtonText="Cancel"
        onConfirm={handleConfirmBulkDelete}
        onCancel={handleCancelBulkDelete}
        isLoading={false}
      />

      <ConfirmationModal
        isOpen={isChatDeleteModalOpen}
        title="Delete Chat"
        message={`Are you sure you want to delete "${chatToDelete?.name}"? This action cannot be undone.`}
        confirmButtonText="Delete"
        cancelButtonText="Cancel"
        onConfirm={handleConfirmChatDelete}
        onCancel={handleCancelChatDelete}
        isLoading={deleteChatSessionPending}
      />
    </>
  );
};
