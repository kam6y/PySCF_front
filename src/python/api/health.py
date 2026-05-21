from fastapi import APIRouter, Request

router = APIRouter()


@router.get('/health')
def health_check(request: Request) -> dict[str, str]:
    version = getattr(request.app.state, 'VERSION', '1.0.0')
    return {
        'status': 'ok',
        'service': 'pyscf-front-api',
        'version': version,
    }
