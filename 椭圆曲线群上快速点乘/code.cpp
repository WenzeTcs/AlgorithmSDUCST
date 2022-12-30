#pragma GCC optimize(2)
#include<bits/stdc++.h>
using namespace std;

#define ll __int128
//map<ll,ll> f;
const ll p=211108170305887;
ll a, b;

#define rep(i,s,t) for(int i=s; i<=t; ++i)

ll qp(ll a, ll x) {
	ll res = 1;
	while(x) {
		if(x&1) res = res*a%p;
		a = a*a%p;
		x >>= 1;
	//	cout << x <<endl;
	}
//	cout <<res<<endl;
	return res;
}int w = 2;
ll inv(ll a) {
	return qp(a, p-2);
}
#define ECPoint pair<ll,ll>
#define xx first
#define yy second

#define x1 A.xx
#define y1 A.yy
#define x2 B.xx
#define y2 B.yy
#define x3 C.xx
#define y3 C.yy
ECPoint add(ECPoint A, ECPoint B) {
	ECPoint C;
	if(x1 == x2 && y1 ==y2) {
		ll k = (3*x1*x1%p + a)*inv(2*y1)%p;
		x3= (k*k%p-2*x1+p+p)%p;
		y3 = (k*((x1-x3+p)%p)%p-y1+p)%p;
	} else {
		ll k = (y2-y1+p)%p * inv((x2-x1+p)%p) %p;
		x3 = (k*k%p - x1 - x2+p+p)%p;
		y3 = (k*((x1-x3+p)%p)%p-y1+p)%p;
	}
	return C;
}
int K[64];
int calc_wNAF(int w, ll k) {
	int i=0;
	while(k) {
		if(k&1) {
			K[i] = (1<<(w-1)) - (k%(1<<w));
			k = k-K[i];
		} else {
			K[i] = 0;
		}
		k >>= 1; ++i;
	}
	return i-1;
}

ECPoint dot_NAF(ECPoint A, ll k) {
	if(k == 1) return A;
	ECPoint Ai[64];
	int i=0;
	Ai[0] = add(A,A);
	Ai[1] = A;
	while(i<(1<<w-1)) {
		++i;
		Ai[i*2+1] = add(Ai[i*2-1], Ai[0]); 
	}
	ECPoint B=A, C=A;
	int m = calc_wNAF(w,k);
	while(--m) {
		C = add(C,C);
		if(K[m]) {
			if(K[m] > 0) C = add(C, Ai[m]);
			else {
				C = add(C,make_pair(Ai[m].xx,(-Ai[m].yy+p)%p));
			}
		}
	}
	return C;
}

ECPoint dot(ECPoint A, ll k) {
	ECPoint B=A, C=A;
	k--;
	while(k) {
		if(k&1) {
			C = add(B, C);
		} 
		B = add(B,B);
		k>>=1;
	}
	return C;
}
map <ECPoint, ll> f;
void print(__int128 x)
{
    if(x<0)
    {
        putchar('-');
        x=-x;
    }
    if(x>9)print(x/10);
    putchar(x%10+'0');
}
int main()
{
	ECPoint P, R;
	P={47815642535808, 116240163507508};
	ll dlog = 0xfffffffffffffffffffffffffffffffffff;
    R={77983503452527, 143728424564583};
	ECPoint C;
	clock_t t0, t1,cnt;
	rep(j,1,10) {
	t0 = clock();
	rep(i,1,10000)
		C = dot(P, dlog);
	t1 = clock(); 
	cnt += t1-t0;
	}
	cout << " Average quick-power-like dot product time cost: "<<((double)cnt/100)/CLOCKS_PER_SEC<<"ms"<<endl;
cnt = 0;
	rep(j,1,10) {
	t0 = clock();
	rep(i,1,10000)
		C = dot_NAF(P, dlog);
	t1 = clock(); 
	cnt += t1-t0;
	}
	cout << "Average wNAF dot product time cost: "<<((double)cnt/100)/CLOCKS_PER_SEC<<"ms"<<endl;
	
//	print(C.xx);
//	cout << "no solution"<<endl;
	return 0;
}
