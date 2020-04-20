functions{
    real maxEigen(matrix M);
    matrix[,] scale_contact(real pWorkOpen, int nAgeGroups){
        matrix[nAgeGroups, nAgeGroups] CM[6,4] = 
        {
            { //BASE
                diag_matrix(rep_vector(1, 16)), //HOME
                diag_matrix(rep_vector(1, 16)), //WORK
                diag_matrix(rep_vector(1, 16)), //SCHOOL
                diag_matrix(rep_vector(1, 16))  //OTHERS
            },
            { //lockdown
                diag_matrix(rep_vector(1, 16)),         //HOME
                diag_matrix(rep_vector(.1, 16)),        //WORK
                diag_matrix(rep_vector(0, 16)),         //SCHOOL
                diag_matrix(rep_vector(.1, 16))         //OTHERS
            },
            { //schoolcloseonly
                diag_matrix(rep_vector(1, 16)),         //HOME
                diag_matrix(rep_vector(.5, 16)),        //WORK
                diag_matrix(rep_vector(1, 16)),         //SCHOOL
                diag_matrix(rep_vector(.1, 16))         //OTHERS
            },
            { //workcloseonly
                diag_matrix(rep_vector(1, 16)),         //HOME
                diag_matrix(rep_vector(.5, 16)),        //WORK
                diag_matrix(rep_vector(1, 16)),         //SCHOOL
                diag_matrix(rep_vector(.1, 16))         //OTHERS
            },
            { //socialdistancing
                diag_matrix(rep_vector(1, 16)),         //HOME
                diag_matrix(rep_vector(pWorkOpen, 16)), //WORK
                diag_matrix(rep_vector(0, 16)),         //SCHOOL
                diag_matrix(rep_vector(.1, 16))         //OTHERS
            },
            { //postoutbreak
                diag_matrix(rep_vector(1, 16)),         //HOME
                diag_matrix(rep_vector(.8, 16)),         //WORK
                diag_matrix(rep_vector(.8, 16)),         //SCHOOL
                diag_matrix(rep_vector(.6, 16))          //OTHERS
            }
        };

        return CM;
    }

    real get_beta(
        real R0,
        int nAgeGroups, 
        real gamma, //removal rate
        int calculate_transmission_probability,
        vector POP,
        matrix[] contact_matrix)
    {
        matrix[nAgeGroups, nAgeGroups] constraints_base[4];
        real TOTALPOP = sum(POP);
        vector[nAgeGroups] pAge = POP/TOTALPOP;

        matrix[nAgeGroups, nAgeGroups] mC[4];
        // vector eigVal;
        real beta;

        matrix[nAgeGroups, nAgeGroups] x;
        matrix[nAgeGroups, nAgeGroups] Csym;
        matrix[nAgeGroups, nAgeGroups] C=rep_matrix(0, nAgeGroups, nAgeGroups);

        for (i in 1:4){
            constraints_base[i] = diag_matrix(rep_vector(1, nAgeGroups));
            x = contact_matrix[i];
            Csym = (x + transpose(x).*((pAge)*transpose(1 ./ pAge)))/2;

            // mC[i] = Csym*constraints_base[i];
            C += Csym*constraints_base[i];
        }

        if (calculate_transmission_probability==1){
            matrix[nAgeGroups, nAgeGroups] M = C;
            for (i in 1:nAgeGroups){
                for (j in 1:nAgeGroups){
                    M[i,j] = pAge[i]/pAge[j]*C[i,j];
                }
            }
            // eigVal =eigenvalues_sym(M);
            // beta = R0*gamma/max(eigenvalues_sym(M));
            beta=R0*gamma/maxEigen(M);
        } else beta = .025;
        
        return beta;
    }

    matrix[] simulate_SEIR(vector POP, vector initialI, //initial I statified by age group 
        real R0, real R0postoutbreak, vector rho, //R0 and R0 post-oubreak
        int nDaySim, real DurInf, real DurLat, matrix[] contact_matrix, 
        int tStartIntenseIntervention, int nIntenseStages, int[] IntenseStagesWeeks, real[] pWorkOpen, 
        int tCloseSchool, int tReopenSchool)
    {
        real gamma = 1.0-exp(-1.0/DurInf);                                       //removal rate
        real alpha = 1.0-exp(-1.0/DurLat);                                       //exposure rate
        int tmax = nDaySim;
        int dt = 1;
        int nSteps = nDaySim;
        int tStopIntenseIntervention = tStartIntenseIntervention + sum(IntenseStagesWeeks)*7;

        int nAgeGroups = dims(POP)[1];
        row_vector[nAgeGroups] numExposed;
        row_vector[nAgeGroups] numInfected;
        row_vector[nAgeGroups] numRecovery;
        row_vector[nAgeGroups] numReported;
        matrix[nSteps, nAgeGroups] S;
        matrix[nSteps, nAgeGroups] I;
        matrix[nSteps, nAgeGroups] E;
        matrix[nSteps, nAgeGroups] R;
        matrix[nSteps, nAgeGroups] H;
        matrix[nSteps, nAgeGroups] lambda;
        matrix[nSteps, nAgeGroups] incidence;
        matrix[nSteps, nAgeGroups] reported;
        matrix[nSteps, nAgeGroups] cumulativeIncidence;
        vector[nSteps] time;
        vector[nSteps] pWork;
        real beta;
        real beta_postfirstwave;
        matrix[nAgeGroups, nAgeGroups] constraintsIntervention[6,4];
        matrix[nAgeGroups, nAgeGroups] CONTRAINTS[4];
        matrix[nAgeGroups, nAgeGroups] C;

        matrix[nSteps, nAgeGroups] output[7];
        
        // Compartments of the models
        S = rep_matrix(0, nSteps, nAgeGroups);
        E = rep_matrix(0, nSteps, nAgeGroups);
        I = rep_matrix(0, nSteps, nAgeGroups);
        R = rep_matrix(0, nSteps, nAgeGroups);
        H = rep_matrix(0, nSteps, nAgeGroups);
        lambda = rep_matrix(0, nSteps,nAgeGroups);
        incidence = rep_matrix(0, nSteps,nAgeGroups);
        reported = rep_matrix(0, nSteps,nAgeGroups);
        cumulativeIncidence = rep_matrix(0, nSteps,nAgeGroups);
        time = rep_vector(0,nSteps);

        I[1,:] = initialI';
        S[1,:] = POP' - I[1,:];

        // Interventions def
        pWork = rep_vector(1, nSteps);
        for (i in 1:dims(IntenseStagesWeeks)[1]){
            int a = i==1 ? tStartIntenseIntervention : tStartIntenseIntervention + sum(IntenseStagesWeeks[1:(i-1)])*7;
            pWork[a:(a+IntenseStagesWeeks[i]*7-1)] = rep_vector(pWorkOpen[i], IntenseStagesWeeks[i]*7);
            // pWork = append_row(pWork, rep_vector(pWorkOpen[i], IntenseStagesWeeks[i]));
        }

        beta = get_beta(R0, nAgeGroups, gamma, 1, POP, contact_matrix);
        beta_postfirstwave = (pWorkOpen[2]<1) ? get_beta(R0postoutbreak, nAgeGroups, gamma, 1, POP, contact_matrix) : beta;
        for (s in 1:(nSteps-1)){
            constraintsIntervention = scale_contact(pWork[s], nAgeGroups);
            
            // Load the intervention scale matrix
            if (time[s] < tCloseSchool || time[s] >= tReopenSchool){
                if (time[s] < tStartIntenseIntervention){
                    CONTRAINTS = constraintsIntervention[1];
                } else if (time[s] >= tStartIntenseIntervention && time[s] < tStopIntenseIntervention){
                    CONTRAINTS = constraintsIntervention[4];
                } else {
                    CONTRAINTS = constraintsIntervention[6];
                }
            } else {
                if (time[s] < tStartIntenseIntervention || time[s] >= tStopIntenseIntervention){
                    CONTRAINTS = constraintsIntervention[2];
                } else {
                    CONTRAINTS = constraintsIntervention[5];
                }     
            }

            C = CONTRAINTS[1]*contact_matrix[1] + 
                CONTRAINTS[2]*contact_matrix[2] + 
                CONTRAINTS[3]*contact_matrix[3] +
                CONTRAINTS[4]*contact_matrix[4];

            // Calculate the force of infection
            lambda[s,:] = time[s] < tStopIntenseIntervention ? (beta*(C*(I[s,:]./POP')'))' : (beta_postfirstwave*(C*(I[s,:]./POP')'))';

            // calculate the number of infections and recoveries between time t and t+dt
            numExposed = lambda[s,:].*S[s,:]*dt; // S to E
            numInfected = alpha*E[s,:]*dt;               // E to I
            numRecovery = gamma*I[s,:]*dt;               // I to R
            numReported = numInfected*rho[s];           // I to H, but not removed from I (remember this is an accumulator function)
    
            S[s+1,:] = S[s,:]-numExposed;
            E[s+1,:] = E[s,:]+numExposed-numInfected;
            I[s+1,:] = I[s,:]+numInfected-numRecovery;
            R[s+1,:] = R[s,:]+numRecovery;
            H[s+1,:] = H[s,:]+numReported;
            incidence[s+1,:] = numInfected/dt;
            reported[s+1,:] = numReported/dt;
            time[s+1] = time[s]+dt;
        }

        output = {S, E, I, R, lambda, incidence, reported};
        return output;
    }
}

data{
    int nDaySim;                                                 //step size
    real mean_DurInf;
    real mean_DurLat;
    real<lower=0> mean_R0; real<lower=0> s_R0;
    real<lower=0> mean_R0postoutbreak; real<lower=0> s_R0postoutbreak;

    int<lower=1> nAgeGroups;
    vector<lower=0>[nAgeGroups] POP;
    vector<lower=0>[nAgeGroups] initialI;
    matrix[nAgeGroups, nAgeGroups] contact_matrix[4];
    vector<lower=0, upper=1>[nDaySim] rho;

    int<lower=0> tStartIntenseIntervention;
    int<lower=0> nIntenseStages;
    int<lower=0> IntenseStageWeeks[nIntenseStages];
    real<lower=0, upper=1> pWorkOpen[nIntenseStages];

    int<lower=0> tCloseSchool;
    int<lower=0> tReopenSchool;
}

parameters{
    real<lower=0, upper=2> lnR0;
    real<upper=2> lnR0postoutbreak;
    // real<lower=0> R0;
    // real<lower=0> R0po
    real<lower=1, upper=30> DurInf; 
    real<lower=1, upper=30> DurLat;
}


model{
    lnR0 ~ normal(log(mean_R0), s_R0)T[0,2];
    lnR0postoutbreak ~ normal(log(mean_R0postoutbreak), s_R0postoutbreak)T[,2];
    DurInf ~ exponential(1.0/mean_DurInf)T[1,]; //Marc said the serial interval is normally distributed?
    DurLat ~ exponential(1.0/mean_DurLat)T[1,];
}

generated quantities{
     matrix[nDaySim, nAgeGroups] SEIR[7] = simulate_SEIR(POP, initialI, exp(lnR0), exp(lnR0postoutbreak), rho,
        nDaySim, DurInf, DurLat, contact_matrix, 
        tStartIntenseIntervention, nIntenseStages, IntenseStageWeeks, pWorkOpen, 
        tCloseSchool, tReopenSchool);
}

