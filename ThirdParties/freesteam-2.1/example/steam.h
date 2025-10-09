typedef struct SteamState_R1_struct{
	double p, T;
} SteamState_R1;

typedef struct SteamState_R2_struct{
	double p, T;
} SteamState_R2;

typedef struct SteamState_R3_struct{
	double rho, T;
} SteamState_R3;

typedef struct SteamState_R4_struct{
	double T, x;
} SteamState_R4;

typedef struct SteamState_struct{
	char region;
	union{
		SteamState_R1 R1;
		SteamState_R2 R2;
		SteamState_R3 R3;
		SteamState_R4 R4;
	};
} SteamState;