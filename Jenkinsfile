pipeline {
    agent { label 'linux' }
    stages {
        stage('Setup') {
            steps {
                echo 'Setting up...'
                echo '- Updating'
                sh 'apt-get update'
                echo '- Installing build tools'
                sh 'apt-get install -y build-essential'
            }
        }
        stage('Build') {
            steps {
                echo 'Building...'
                sh 'make -C experiments/xor/'
            }
        }
        stage('Test') {
            steps {
                echo 'Testing...'
                sh './experiments/xor/xor'
                sh './experiments/xor/xor_independent'
            }
        }
    }
}