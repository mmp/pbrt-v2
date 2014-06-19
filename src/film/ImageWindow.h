#pragma once
#include <assert.h>
class ImageWindow
{
	HWND m_hWnd;
	int m_Width;
	int m_Height;
	HANDLE	m_hThread;
	DWORD * m_Bits;
	static DWORD WINAPI ThreadProc(LPVOID lpParameter);
public:
	ImageWindow(int width, int height);
	~ImageWindow();
	__inline void SetPixel(int x, int y,DWORD c){
		assert(x >= 0 && x < m_Width);
		assert(y >= 0 && y < m_Height);
		m_Bits[x + y * m_Width] = c;
	}
	void UpdateWindow();
};