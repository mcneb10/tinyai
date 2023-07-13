pipeline {
    agent linux
    stages {
        stage('Setup') {
            steps {
                echo 'Updating...'
                sh 'apt update'
                echo '- Installing build tools'
                sh 'apt install build-essential'
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