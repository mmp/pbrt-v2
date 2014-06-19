#include "stdafx.h"
#include "ImageWindow.h"

#include <tchar.h>
static TCHAR szWindowClass[1024] = { _T("SimpleTestWindowClass") };
static TCHAR szTitle[1024] = { _T("Window Title") };
static DWORD dwStyle = WS_OVERLAPPED|WS_CAPTION|WS_MINIMIZEBOX|WS_SYSMENU;
static DWORD dwStyleEx = 0;
static LRESULT CALLBACK WndProc(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam);
static
ATOM MyRegisterClass()
{
	HINSTANCE hInstance = GetModuleHandle(NULL);
	WNDCLASSEX wcex;
	
	wcex.cbSize = sizeof(WNDCLASSEX);

	wcex.style = CS_HREDRAW | CS_VREDRAW;
	wcex.lpfnWndProc = WndProc;
	wcex.cbClsExtra = 0;
	wcex.cbWndExtra = 0;
	wcex.hInstance = hInstance;
	wcex.hIcon = 0;// LoadIcon(hInstance, MAKEINTRESOURCE(IDI_SIMPLEWINDOW));
	wcex.hCursor = LoadCursor(NULL, IDC_ARROW);
	wcex.hbrBackground = (HBRUSH)(COLOR_WINDOW + 1);
	wcex.lpszMenuName = NULL;//MAKEINTRESOURCE(IDC_SIMPLEWINDOW);
	wcex.lpszClassName = szWindowClass;
	wcex.hIconSm = NULL;//LoadIcon(wcex.hInstance, MAKEINTRESOURCE(IDI_SMALL));

	return RegisterClassEx(&wcex);
}
static
HWND InitInstance()
{
	HINSTANCE hInstance = GetModuleHandle(NULL);
	HWND hWnd;

	//hInst = hInstance; // Store instance handle in our global variable

	hWnd = CreateWindow(szWindowClass, szTitle, dwStyle,
		CW_USEDEFAULT, 0, CW_USEDEFAULT, 0, NULL, NULL, hInstance, NULL);

	if (!hWnd)
	{
		return NULL;
	}

	ShowWindow(hWnd, SW_SHOW);
	UpdateWindow(hWnd);

	return hWnd;
}

#define WM_UPDATEBITMAP (WM_USER + 100)
struct UpdateBitmapInfo{
	int width;
	int height;
	DWORD * bits;
};
static
LRESULT CALLBACK WndProc(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam)
{
	PAINTSTRUCT ps;
	HDC hdc;
	HBITMAP hBitmap = NULL;
	UpdateBitmapInfo * ubInfo = NULL;

	switch (message)
	{
	case WM_UPDATEBITMAP:
		hdc = GetDC(hWnd);
		ubInfo = (UpdateBitmapInfo*)wParam;
		hBitmap = (HBITMAP)GetProp(hWnd, _T("my_bitmap"));
		if (hBitmap == NULL)
		{
			hBitmap = CreateCompatibleBitmap(hdc, ubInfo->width, ubInfo->height);
			SetProp(hWnd, _T("my_bitmap"), hBitmap);
		}
		{
			BITMAPINFO bmInfo = { 0 };
			bmInfo.bmiHeader.biSize = sizeof(bmInfo.bmiHeader);
			bmInfo.bmiHeader.biWidth = ubInfo->width;
			bmInfo.bmiHeader.biHeight = -ubInfo->height;
			bmInfo.bmiHeader.biPlanes = 1;
			bmInfo.bmiHeader.biBitCount = 32;
			bmInfo.bmiHeader.biCompression = BI_RGB;
			bmInfo.bmiHeader.biSizeImage = 0;
			SetDIBits(hdc, hBitmap, 0, ubInfo->height, ubInfo->bits, &bmInfo, 0);
		}
		ReleaseDC(hWnd,hdc);
		InvalidateRect(hWnd, NULL, FALSE);
		break;
	case WM_PAINT:
		hdc = BeginPaint(hWnd, &ps);
		hBitmap = (HBITMAP)GetProp(hWnd, _T("my_bitmap"));
		if (hBitmap != NULL)
		{
			BITMAP bm;
			GetObject(hBitmap, sizeof(bm), &bm);
			HDC hmem = CreateCompatibleDC(hdc);
			HBITMAP old = (HBITMAP)SelectObject(hmem, hBitmap);
			BitBlt(hdc, 0, 0, bm.bmWidth,bm.bmHeight, hmem, 0, 0, SRCCOPY);
			SelectObject(hmem, old);
			DeleteDC(hmem);
		}
		EndPaint(hWnd, &ps);
		break;
	case WM_DESTROY:
		hBitmap = (HBITMAP)GetProp(hWnd, _T("my_bitmap"));
		if (hBitmap)
		{
			DeleteObject(hBitmap);
			hBitmap = NULL;
			RemoveProp(hWnd, _T("my_bitmap"));

		}
		PostQuitMessage(0);
		break;
	default:
		return DefWindowProc(hWnd, message, wParam, lParam);
	}
	return 0;
}

DWORD WINAPI ImageWindow::ThreadProc(LPVOID lpParameter)
{
	ImageWindow * self = (ImageWindow*)lpParameter;
	MSG msg;
	
	MyRegisterClass();
	HWND hWnd = self->m_hWnd = InitInstance();
	if (!hWnd)
	{
		return 0;
	}

	{
		RECT rc;
		GetWindowRect(hWnd, &rc);
		rc.right = rc.left + self->m_Width;
		rc.bottom = rc.top + self->m_Height;
		AdjustWindowRect(&rc, dwStyle, FALSE);
		SetWindowPos(hWnd, NULL, 0, 0, rc.right - rc.left, rc.bottom - rc.top, SWP_NOMOVE | SWP_NOZORDER);
	}

	while (GetMessage(&msg, NULL, 0, 0))
	{
		TranslateMessage(&msg);
		DispatchMessage(&msg);
	}
	self->m_hWnd = NULL;
	return 1;
}

ImageWindow::ImageWindow(int width, int height)
{
	m_hWnd = NULL;
	m_Width = width;
	m_Height = height;
	m_Bits = (DWORD*)malloc(width * height * sizeof(DWORD));
	memset(m_Bits,0,width * height * sizeof(DWORD));
	DWORD id;
	m_hThread = CreateThread(NULL, NULL, &ThreadProc, this, NULL, &id);
	
	DWORD ret;
	while (m_hWnd == NULL && GetExitCodeThread(m_hThread, &ret) && ret == STILL_ACTIVE)
	{
		Sleep(10);
	}
	UpdateWindow();
}
ImageWindow::~ImageWindow()
{
	WaitForSingleObject(m_hThread, INFINITE);
	free(m_Bits);
	CloseHandle(m_hThread);
}

void ImageWindow::UpdateWindow()
{
	UpdateBitmapInfo ub;
	ub.width = m_Width;
	ub.height = m_Height;
	ub.bits = m_Bits;
	DWORD ret;
	if (m_hWnd && GetExitCodeThread(m_hThread, &ret) && ret == STILL_ACTIVE)
	{
		SendMessage(m_hWnd, WM_UPDATEBITMAP, (WPARAM)&ub, 0);
	}
	
}